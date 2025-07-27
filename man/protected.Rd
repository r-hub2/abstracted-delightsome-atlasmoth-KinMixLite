\name{protected}
\alias{protected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Catch numerical errors, and return -Inf}
\description{Attempts to catch numerical erros in evaluating the expression \code{x}, delivering a default result instead of NaN or other failures}
\usage{
protected(x,default=-Inf) 
}
\arguments{
  \item{x}{expression to be evaluated, typically the log-likelihood of a modified mixture model}
  \item{default}{value to be delivered if numerical errors are encountered}
}
\value{Returns -Inf in case of error, otherwise the value of \code{x}}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
