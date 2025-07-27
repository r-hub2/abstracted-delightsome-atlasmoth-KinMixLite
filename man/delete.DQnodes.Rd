\name{delete.DQnodes}
\alias{delete.DQnodes}
\title{Delete D and Q dummy nodes from all Bayes nets in mixture}
\description{Delete D and Q dummy nodes and associated edges from all Bayes nets in mixture, to save space; these nodes would 
only be needed for specific follow-up analyses}
\usage{
delete.DQnodes(mixture,which="DQ")
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{which}{character string}
}
\details{The function removes the D and/or Q nodes from the DNAmixture object, depending on whether \code{which} includes "D", "Q" or both}
\value{No value is returned, the function is called for its side effect}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

\examples{
data(test2data)
data(NGMDyes)

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,dyes=list(NGMDyes),
	triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
replace.tables.for.UAF(mixD,40)
size(mixD)
}
