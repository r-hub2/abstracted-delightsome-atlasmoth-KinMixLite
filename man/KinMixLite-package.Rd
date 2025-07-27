\name{KinMixLite-package}
\alias{KinMixLite-package}
\alias{KinMixLite}
\docType{package}
\title{
\packageTitle{KinMixLite}
}
\description{
\packageDescription{KinMixLite}
}
\details{
This package is a toolkit for inference about mixtures and familial relationships, either between contributors or between a contributor and other typed individuals. It extends the functionality proposed in Green and Mortera (2017) by allowing more general relationships, specified in general by an IBD pattern distribution - the generalisation to more than two individuals of the coefficients of identity of Jacquard (1974). Details are in the paper by Green and Mortera (2021). \code{KinMixLite} extends the capability of the \pkg{DNAmixturesLite} package, and intimately relies on that package; as with that package, instead of the \pkg{RHugin} package, it uses \pkg{gRaven} and \pkg{gRain} for Bayes Net calculations. This version implements the ALN, MBN and WLR as well as RPT methods; see Green and Mortera (2017). 
}
\author{
Maintainer: \packageMaintainer{KinMixLite}
}
\references{
Green, P. J. and Mortera, J. (2017). Paternity testing and other inference about relationships from DNA mixtures.
\emph{Forensic Science International: Genetics}. <doi:10.1016/j.fsigen.2017.02.001>. 

Green, P. J. and Mortera, J. (2021). Inference about complex relationships using peak height data from DNA mixtures.
\emph{Applied Statistics}. <doi:10.1111/rssc.12498>.

Jacquard, A. (1974) 
\emph{The genetic structure of populations}. Springer-Verlag.
}
\keyword{ package }
\seealso{
\code{\link[DNAmixturesLite:DNAmixturesLite-package]{DNAmixturesLite}}
}
\section{Formats}{
See \code{\link{formats}} for formats of the various data objects used in this package.
}

\examples{
require(ribd)
data(test2data)
data(NGMDyes)

C<-50

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,'parent',list(c=Cgt),targets=c('f','c'),contrib='f') 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture, where contributors are siblings

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.6,U2=0.3,U3=0.1)))
baseline<-logL(mixD)(pars)

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,'sibs',targets=c('b1','b2'),contribs=c('b1','b2')) 
log10LR<-(protected(logL(mixD)(pars))-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
