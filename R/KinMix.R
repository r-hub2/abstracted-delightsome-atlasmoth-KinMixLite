KinMix<-function(data,k,C,database,K=character(0),reference.profiles = NULL, 
	contribs=NULL,typed.gts=NULL,IBD=NULL,targets=NULL,pars=NULL,mle=FALSE,
	dir=character(0),domainnamelist=NULL,
	load=FALSE,write=FALSE,dyes=NULL, 
	triangulate=TRUE,compile=TRUE,compress=TRUE,use.order=TRUE)
{
checkdata(data,database,typed.gts=NULL,reference.profiles=NULL)

typed<-unique(names(typed.gts))
if(length(typed)!=length(typed.gts)) stop('names of typed.gts must be distinct')
targets<-unique(targets)
IBD<-as.IBD(IBD)
if(length(targets)*2!=ncol(IBD$patt)) stop('targets and IBD inconsistent')

z<-KMmodel(k,typed,targets,contribs,rp=TRUE,verbose=TRUE)


xprofiles<-z$move
# merge both sorts of known contributors
if(!is.null(reference.profiles)){
	if(length(xprofiles)>0) {
		extra<-make.profiles(typed.gts[xprofiles])
		reference.profiles<-merge(reference.profiles,extra,all=TRUE)
		reference.profiles[is.na(reference.profiles)]<-0
	}
} else {
	if(length(xprofiles)>0) {
		reference.profiles<-make.profiles(typed.gts[xprofiles])
	}
}
xprofiles<-c(K,xprofiles)
xcontribs<-z$contribs
xtyped<-z$typed
if(is.null(z$typed)) typed.gts<-NULL else typed.gts<-typed.gts[z$typed]

xignore<-setdiff(typed,union(contribs,targets))

#cat('[[ diagnostic:\nprofiles:',xprofiles,'\n')
#cat('contribs:',xcontribs,'\n')
#cat('typed:   ',xtyped,'\n')
#cat('targets: ',targets,'\n')
#cat('ignore:  ',xignore,'\n]]\n')


if(0==length(xcontribs)){
if(0==length(xprofiles)) {

#cat('** DNAmixture, no profiles **\n')
mix<-DNAmixture(data,k,C=C,database=database)

}else{

#cat('** DNAmixture, with profiles **\n')
mix<-DNAmixture(data,k,K=xprofiles,reference.profiles=reference.profiles,C=C,database=database)

}
}else{
if(0==length(xprofiles)&&is.null(reference.profiles)) {

#cat('** DNAmixture with rpt.IBD, no profiles **\n')
mix<-DNAmixture(data,k,C=C,database=database,triangulate=FALSE)
rpt.IBD(mix,IBD,typed.gts=typed.gts[xtyped],targets=targets,contribs=xcontribs,quiet=TRUE) 

}else{

#cat('** DNAmixture, with profiles and rpt.IBD **\n')
mix<-DNAmixture(data,k,K=xprofiles,reference.profiles=reference.profiles,C=C,database=database,triangulate=FALSE)
rpt.IBD(mix,IBD,typed.gts=typed.gts[xtyped],targets=targets,contribs=xcontribs,quiet=TRUE) 
}
}

# following assumes only 1 trace!

if(!is.null(pars)) 
{
# nx<-c(xprofiles,xcontribs); if(length(nx)<k) nx<-c(nx,paste0("U",(1:(k-length(nx)))))
new<-c(xprofiles); if(length(new)<k) new<-c(new,paste0("U",(1:(k-length(new)))))
parsi<-pars
now<-names(parsi[[4]])
if(is.null(now)) {
# phi unnamed case
	names(parsi[[4]])<-new 
	cat('phi named in pars: ',paste(names(parsi[[4]]), collape=' '),'\n')
	} else {
# phi may need to be renamed case
	if(!all(now%in%new))
	{
	names(parsi[[4]])[!(now%in%new)]<-setdiff(new,now)
	cat('phi renamed in pars: ',paste(now, collapse=' '),' to ',paste(names(parsi[[4]]), collape=' '),'\n')
	}
}
#cat('[[ diagnostic\ninternal naming:\n')
#print(parsi)
#cat(']]\n')

if(mle) {
	mlm<-mixML(mix,parsi)
	mix$pars<-mix$mle<-mlm$mle
	mix$logL<-mlm$lik
} else {
	mix$logL<-logL(mix)(parsi)
	mix$pars<-parsi
}
}
mix
}

make.profiles<-function(typed.gts)
{
nt<-names(typed.gts)
profiles<-make.profile(typed.gts[[1]],nt[1])
if(length(typed.gts)>1) for(i in 2:length(typed.gts)) 
	profiles<-merge(profiles,make.profile(typed.gts[[i]],nt[i]),all=TRUE)
profiles[is.na(profiles)]<-0
profiles
}
