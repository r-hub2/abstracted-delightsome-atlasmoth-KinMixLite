KMmodel<-function(k,typed,targets,contribs,rp=FALSE,verbose=FALSE)
{
typed<-unique(typed)
targets<-unique(targets)
contribs<-unique(contribs)
xcontribs<-contribs
#
dtarg<-setdiff(targets,union(typed,contribs))
if(verbose) if(length(dtarg)>0) cat(dtarg,'dropped from targets (no effect)\n')
xtargets<-setdiff(targets,dtarg)
#
dtyped<-setdiff(typed,union(targets,contribs))
if(verbose) if(length(dtyped)>0) cat(dtyped,'dropped from typed (no effect)\n')
xtyped<-setdiff(typed,dtyped)
#
#if(length(contribs)>k)
#{
#dcont<-contribs[-(1:k)]
#if(verbose) cat(dcont,'dropped from contribs, given value of k\n')
#xcontribs<-contribs[1:k]
#} else if(length(contribs)<k) 
#{
#xcontribs<-c(contribs,paste0('U',1:(k-length(contribs))))
#if(verbose) cat('given value of k, contribs extended to',xcontribs,'\n')} else xcontribs<-contribs
#
move<-setdiff(intersect(typed,contribs),targets)
if(length(move)>0) {
	if(rp) {
	xtyped<-setdiff(xtyped,move)
	xcontribs<-setdiff(xcontribs,move)
	cat(move,'will be added to ref prof\n',
	'having been dropped from typed and contribs\n')
	} else {
	xtargets<-union(xtargets,move)
	cat(move,'added to targets and should be added to IBD\n')
	}
}
#
cat('Model:\n      typed:',xtyped,'\n    targets:',xtargets,'\n   contribs:',xcontribs,'\n       move:',move,'\n')
invisible(list(typed=xtyped,targets=xtargets,contribs=xcontribs,move=move))
}
