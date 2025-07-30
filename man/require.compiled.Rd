\name{require.compiled}
\alias{require.compiled}
\title{Force compilation of all BNs in a DNA mixture model}
\description{Scan all Bayes nets in mixture, and compile any that are not already compiled}
\usage{
require.compiled(mixture)
}
\arguments{
  \item{mixture}{A DNAmixture object}
}
\value{No value is returned, the function is called for its side effect}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}


\examples{
data(test2data)
data(NGMDyes)

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,dyes=list(NGMDyes),
	triangulate=FALSE,compile=FALSE)
replace.tables.for.UAF(mixD,40,compile=FALSE)
require.compiled(mixD)
}
