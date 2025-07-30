\name{add.motherchild.likd.node}
\alias{add.motherchild.likd.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replace CPTs for mixture contributor a Father, given Child and Mother genotypes, by ALN method}
\description{loop over markers, and alleles within markers
to create node Rlikd for relative likelihood for individual i,
for paternity model with Mother and Child genotyped
then compile all domains. Implements method ALN.
}
\usage{
add.motherchild.likd.node(mixture,Cgt,Mgt,db,ind=1)
}
\arguments{
  \item{mixture}{A DNAmixture object}
  \item{Cgt}{Child's genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
  \item{Mgt}{Mother's genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
  \item{db}{Allele frequency database}
  \item{ind}{Index of contributor regarded as Father: which `unknown' contributor are we modelling by amending his CPTs?}
}
\value{No value is returned, the function is called for its side effect}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

mixD3<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
cgtcaca<-gt2aca(mixD3,Cgt) 
add.motherchild.likd.node(mixD3,Cgt,Mgt,db,1)
log10LR<-(logLX(mixD3,
	expr.make.findings(list(
	list('Male',ind=1),
	list('Rlikd',aca='cgtcaca',cgt='Cgt',evid='Revid')
	))
	)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n') 
}
