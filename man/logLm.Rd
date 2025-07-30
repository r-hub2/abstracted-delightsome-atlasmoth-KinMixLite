\name{logLm}
\alias{logLm}
\title{Log likelihood for mixture model, by marker}
\description{The function logLm is used to produce a log likelihood function for one or more traces of DNA, delivering results on a per-marker basis.}
\usage{
logLm(mixture, presence.only = FALSE, initialize = TRUE)
}
\arguments{
\item{mixture}{A \code{\link[DNAmixturesLite]{DNAmixture}} model.}
\item{presence.only}{Set to TRUE to ignore peak heights and
evaluate the likelihood function considering peak presence and
absence (heights above and below threshold) only. Defaults to
FALSE.}
\item{initialize}{By default all entered
evidence is removed from the networks in \code{object}. Setting
\code{initialize = FALSE} should be done with care, and it is up
to the user to ensure that the likelihood being computed is meaningful.}
}
\value{
A function, which takes a \code{\link[DNAmixturesLite]{mixpar}} model parameter as argument, and delivers a vector of per-marker loglikelihood values.
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\seealso{
See also \code{\link[DNAmixturesLite]{logL}}.
}
\examples{
data(test2data)
data(NGMDyes)

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logLm(mixD)(pars)

## Fit 2-person mixture model in which contributor 1 is father of a typed individual Cgt 
## with mother Mgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,,list(Mgt,Cgt)) 
log10LR<-(logLm(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,' sum:',sum(log10LR),'\n')
}
