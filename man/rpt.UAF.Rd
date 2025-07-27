\name{rpt.UAF}
\alias{rpt.UAF}
\alias{replace.tables.for.UAF}
\title{Replace CPTs in a DNA mixture to model uncertain allele frequencies}
\description{
Replace CPTs in a DNA mixture to model uncertainty in allele frequencies
}
\usage{
replace.tables.for.UAF(mixture, M, compile = TRUE)

rpt.UAF(mixture, M, compile = TRUE)
}
\arguments{
  \item{mixture}{DNAmixtures object created by previous call to \code{DNAmixture} with \cr
\code{triangulate=FALSE,compile=FALSE}}
  \item{M}{Size of allele frequency database}
  \item{compile}{logical flag: should mixture object be compiled on exit?}
}
\value{No value is returned, the function is called for its side effect}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\examples{
data(test2data)
data(NGMDyes)

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

## Fit 2-person mixture model under assumption that database size was only 40

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
replace.tables.for.UAF(mixD,40)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
