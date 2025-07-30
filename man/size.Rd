\name{size}
\alias{size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate and display total size of BN tables for a DNA mixture}
\description{Calculate and display total size of BN tables for a DNA mixture}
\usage{
size(mixture)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
}
\value{Returns total size, typically to be printed by bespoke method}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)
data(NGMDyes)

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
size(mixD)
}
