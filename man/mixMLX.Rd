\name{mixMLX}
\alias{mixMLX}
\title{Maximises mixture likelihood when extra findings present}
\description{Replacement for mixML in DNAmixturesLite that calls logLX instead of logL.}
\usage{
mixMLX(mixture, expr.extra.findings, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, initialize = FALSE, ...)
}
\arguments{
  \item{mixture}{Compiled DNAmixture object.}
  \item{expr.extra.findings}{expression containing the extra findings}
  \item{pars}{Parameters, in \code{mixpar} format.}
  \item{constraints}{as in \code{mixML}}
  \item{phi.eq}{as in \code{mixML}}
  \item{val}{as in \code{mixML}}
  \item{trace}{as in \code{mixML}}
  \item{order.unknowns}{as in \code{mixML}}
  \item{initialize}{should all entered evidence be removed from the networks in \code{mixture}}
  \item{...}{as in \code{mixML}}
}
\value{A list containing
\describe{
\item{mle}{The (suggested) MLE.}
\item{lik}{The log of the likelihood (log e).}
}
as well as the output from the optimisation.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link[DNAmixturesLite]{mixML}}.
}
\examples{
data(test2data)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db)

# find MLE's and maximised likelihood
# adding evidence individual 1 is Male 

expr.extra.findings<-expr.make.findings(list(list('Male',ind=1)))

startpar<-mixpar(rho=list(60),eta=list(24),xi=list(0.16),phi=list(c(U1=0.75,U2=0.25)))
mlDM<-mixMLX(mixD,expr.extra.findings,startpar,trace=FALSE) 
pars<-mlDM$mle
cat('\nBaseline model maximised likelihood:',mlDM$lik,'\n')
cat('and MLEs:\n')
print(mlDM$mle)
}
