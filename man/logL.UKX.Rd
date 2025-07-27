\name{logL.UKX}
\alias{logL.UKX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates mixture log likelihood for unknown contributors with extra findings}
\description{Replacement for logL.UK in DNAmixtures that calls \code{extra.findings} immediately before 
propagating all findings and returning the normalising constant for the network.}
\usage{logL.UKX(mixture, expr.extra.findings, initialize = FALSE)}
\arguments{
  \item{mixture}{Compiled DNAmixture object.}
  \item{expr.extra.findings}{expression containing the extra findings}
  \item{initialize}{should all entered evidence be removed from the networks in \code{mixture}}
}
\value{The log likelihood.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}


\seealso{
See also \code{\link[DNAmixturesLite]{logL.UK}}.
}
\examples{
data(test2data)

# set threshold C
C<-0.001

pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.9,U2=0.1)))

mixMBN<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
cgtcaca<-gt2aca(mixMBN,Cgt)
add.child.meiosis.nodes(mixMBN,cgtcaca,1)
logL.UKX(mixMBN,
	expr.make.findings(list(
	list('Male',ind=1),
	list('Caca',aca='cgtcaca')
	)))(pars)

}
