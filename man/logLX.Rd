\name{logLX}
\alias{logLX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates mixture log likelihood when extra findings present}
\description{Replacement for logL in DNAmixturesLite that calls calls LogL.UKX instead of logL.UK.}
\usage{logLX(mixture, expr.extra.findings, presence.only = FALSE, initialize = FALSE)}
\arguments{
  \item{mixture}{Compiled DNAmixture object.}
  \item{expr.extra.findings}{expression containing the extra findings}
  \item{presence.only}{Set to TRUE to ignore peak heights and evaluate the likelihood function considering peak presence and absence (heights above and below threshold) only. Defaults to FALSE}
  \item{initialize}{should all entered evidence be removed from the networks in \code{mixture}}
}
\value{The log likelihood.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\seealso{
See also \code{\link[DNAmixturesLite]{logL}}.
}
\examples{
data(test2data)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

mixMBN<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
cgtcaca<-gt2aca(mixMBN,Cgt)
add.child.meiosis.nodes(mixMBN,cgtcaca,1)
log10LR<-(logLX(mixMBN,
	expr.make.findings(list(
	list('Male',ind=1),
	list('Caca',aca='cgtcaca')
	))
	)(pars)-attr(cgtcaca,'logGt')-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
