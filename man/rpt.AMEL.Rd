\name{rpt.AMEL}
\alias{rpt.AMEL}
\title{Replace CPTs for AMEL marker in a DNA mixture}
\description{Used after a call to \code{DNAmixture} with \code{compile=FALSE,triangulate=FALSE}, this function replaces the CPTs for the genotype allele count arrays for the AMEL marker in a DNA mixture to specify sex of contributors}
\usage{
rpt.AMEL(mixture,sex,compile=TRUE)
}
\arguments{
  \item{mixture}{A DNAmixture object}
  \item{sex}{Integer vector, sex of each contributor}
  \item{compile}{Logical, should BN be compiled after modification?}
}

\details{
The sex of each contributor is coded as in \code{pedtools}, namely 0=unspecified, 1=male, 2=female.
}
\value{No value is returned, the function is called for its side effect}

\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)
data(NGMDyes)

mixD<-DNAmixture(list(epg),k=3,C=list(0.001),database=db,dyes=list(NGMDyes),
        triangulate=FALSE,compile=FALSE)
rpt.AMEL(mixD,c(1,2,0))  # the 3 contributors are male, female, and unspecified,respectively.
}
