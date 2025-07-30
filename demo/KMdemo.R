data(test2data)

t1<-Sys.time()

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)
baseline

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
mlD <- mixML(mixD, pars)
print(mlD$mle)
pars<-mlD$mle

baseline<-logL(mixD)(pars)
baseline

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt, making maximum use
## of defaults

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,,list(Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, more explicitly

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,contribs=c('F'),typed.gts=list(C=Cgt),IBD='parent-child',targets=c('F','C')) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, using wrapper function KinMix

mixD<-KinMix(list(epg),k=2,C=list(0.001),database=db,
	contribs=c('F'),typed.gts=list(C=Cgt),IBD='parent-child',targets=c('F','C'),
	pars=pars)
log10LR<-(mixD$logL-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributor 1 is father of a typed individual Cgt 
## with mother Mgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,IBD='trio',list(Mgt,Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, using wrapper function KinMix

mixD<-KinMix(list(epg),k=2,C=list(0.001),database=db,
	contribs=c('F'),typed.gts=list(M=Mgt,C=Cgt),IBD='trio',targets=c('F','M','C'),
	pars=pars)
log10LR<-(mixD$logL-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributors are two parents of a child with 
## genotype Cgt, and a parent of one of them has genotype Rgt. Note the manual encoding of allele 
## labels to reduce the complexity of the IBD pattern distribution IBD.

IBD<-list(patt=rbind(c(1,3,2,4,1,2,1,5),c(1,3,2,4,1,2,3,5)))

mixD1<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD1,contribs=c('M','F'),typed.gts=list(C=Cgt,GM=Rgt),
	IBD=IBD,targets=c('M','F','C','GM')) 
log10LR<-(logL(mixD1)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, using KinMix, and constructing IBD using as.IBD applied to a matrix of trios defining the pedigree

IBD<-as.IBD(matrix(c('M','GM','GF','C','M','F'),2,3,byrow=TRUE),targets=c('GM','M','F','C'))
mixD2<-KinMix(list(epg),k=2,C=list(0.001),database=db,
	contribs=c('M','F'),typed.gts=list(C=Cgt,GM=Rgt),
	IBD=IBD,targets=c('GM','M','F','C'),pars=pars)
log10LR<-(mixD2$logL-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## show 'deconvolution'

summary(map.genotypes(mixD2,pmin=0.001,U=1:2))



## fit 3-contributor mixture, including one fixed contributor R with genotype profile Rgt

mixD<-DNAmixture(list(epg),k=3,K='R',reference.profiles=make.profile(Rgt,'R'),C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.45,U2=0.28,R=0.27)))
mlD <- mixML(mixD, pars)
print(mlD$mle)
pars<-mlD$mle
print(logL(mixD)(pars))

## now with putative father as one of the 2 unknown contributors

mixD<-DNAmixture(list(epg),k=3,K='R',reference.profiles=make.profile(Rgt,'R'),C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,contribs=c('F'),typed.gts=list(C=Cgt),IBD='parent-child',targets=c('F','C')) 
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.45,U2=0.28,R=0.27)))
mlD <- mixML(mixD, pars)
print(mlD$mle)
pars<-mlD$mle
print(logL(mixD)(pars))

## the same using wrapper function KinMix

mixD<-KinMix(list(epg),k=3,K='R',reference.profiles=make.profile(Rgt,'R'),C=list(0.001),database=db,
	contribs=c('F'),typed.gts=list(C=Cgt),IBD='parent-child',targets=c('F','C'),pars=pars,mle=TRUE)
print(mixD$mle)
print(mixD$logL)

## the same using wrapper function KinMix, including R among the typed genotypes instead

mixD<-KinMix(list(epg),k=3,C=list(0.001),database=db,
	contribs=c('F','R'),typed.gts=list(C=Cgt,R=Rgt),IBD='parent-child',targets=c('F','C'),pars=pars,mle=TRUE)
print(mixD$mle)
print(mixD$logL)


