\name{protected.mixML}
\alias{protected.mixML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Protect against numerical errors in maximum likelihood computation}
\description{Attempts to catch numerical errors in maximum likelihood computation, by replacing logL values by a default value instead of NaN or other failures}
\usage{
protected.mixML(mixture, pars, constraints = NULL, phi.eq = FALSE, val = NULL,
  trace = FALSE, order.unknowns = TRUE, default=-999999, ...)
}
\arguments{
  \item{mixture}{A DNAmixture object.}
  \item{pars}{A mixpar parameter used as a starting value for the optimisation.}
  \item{constraints}{Equality constraint function as function of an array of parameters.}
  \item{phi.eq}{Should the mixture proportions be the same for all traces? Defaults to FALSE.}
  \item{val}{Vector of values to be satisfied for the equality constraints.}
  \item{trace}{Print the evaluations of the likelihood-function during optimisation?}
  \item{order.unknowns}{Should unknown contributors be ordered according to decreasing contributions? Defaults to TRUE.}
  \item{...}{Further arguments to be passed on to solnp.}
  \item{default}{value of logL to be used if numerical errors are encountered}
}
\value{
A list containing
   \item{mle}{The (suggested) MLE.}
   \item{lik}{The log of the likelihood (log e).}
as well as the output from the optimisation.
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
