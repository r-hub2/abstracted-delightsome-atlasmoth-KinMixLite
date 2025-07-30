\name{add.relative.likd.node}
\alias{add.relative.likd.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replace CPTs for mixture contributor a Father, given Child genotype, by ALN method}
\description{loop over markers, and alleles within markers
to create node Rlikd for relative likelihood for individual i,
for paternity model with only Child genotyped
then compile all domains. Implements method ALN.
}
\usage{
add.relative.likd.node(mixture,aca,ind=1)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{aca}{Child's genotype profile as an allele count array}
  \item{ind}{Index of contributor regarded as Parent (or Child): which `unknown' contributor are we modelling by amending his/her CPTs?}
}
\value{No value is returned, the function is called for its side effect}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

mixALN<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
cgtcaca<-gt2aca(mixALN,Cgt) 
add.relative.likd.node(mixALN,cgtcaca,1)
log10LR<-(logLX(mixALN,
	expr.make.findings(list(
	list('Male',ind=1),
	list('Rlikd',aca='cgtcaca',cgt='Cgt',evid='Revid')
	))
	)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n') 
}
