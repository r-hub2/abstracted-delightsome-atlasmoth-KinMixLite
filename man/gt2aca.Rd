\name{gt2aca}
\alias{gt2aca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converts genotype profiles to allele count arrays}
\description{Returns list of vectors of allele counts corresponding to genotype profile in gt}
\usage{
gt2aca(mixture,gt,eps=0)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{gt}{Genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
  \item{eps}{If non-zero, the function creates the output allele count arrays in a different format, that mitigates subsequent propagation errors in some situations. Instead of a vector of allele counts, each element of the list is a matrix with 3 columns, corresponding to allele counts 0, 1 and 2, with entries 1 or \code{eps}.}
}
\value{Returns list of vectors of allele counts. The log probability for the genotype is returned in its attribute '\code{logGt}'.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)
data(NGMDyes)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes))
cgtcaca<-gt2aca(mixD,Cgt)

print(Cgt)
print(cgtcaca)
}
