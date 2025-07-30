\name{convertIBD}
\alias{as.IBD}
\alias{convertIBD}
\title{Convert relationship information to IBD pattern distribution}
\description{Construct IBD pattern distribution from one of several alternative representations
of multi-person condensed coefficients of identity}
\usage{
as.IBD(x='sibs', targets=NULL, ped=FALSE) 
convertIBD(x='sibs', targets=NULL, ped=FALSE)
}
\arguments{
  \item{x}{A string, a vector of length 3 or 9, a list with components \code{pr} and \code{patt}, or a list with two components, a pedigree and a vector of target id's; see Details}
  \item{targets}{character vector of individual tags}
  \item{ped}{logical, should complete pedigree be added as an attribute to the output, if available?}
}

\details{
Possible formats for the input \code{x} are:
\enumerate{
\item certain verbal mnemonics; currently one of the following (or an unambiguous partial match):
c('sibs','parent-child','half-sibs',
	'cousins','half-cousins','second-cousins',
	'double-first-cousins',
	'quadruple-half-first-cousins',
	'3cousins-cyclic','3cousins-star','trio')
\item a vector of 3 kappas
\item a vector of 9 Deltas
\item a list with matrix or vector valued component \code{patt}, with or without component \code{pr}
\item a list with 2 components, the first being a pedigree in the sense of the \code{pedtools} package, the second a vector of target id's
\item a 3-column character matrix of individual tags, each row denoting a child/mother/father triple - an alternative compact representation of a pedigree 
}
}

\value{IBD pattern distribution as a list with components \code{pr} and \code{patt}}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}



\seealso{
\code{\link[pedtools]{pedtools}}, \code{\link{formats}}


}
\examples{
data(test2data)

IBD<-convertIBD('parent-child')

IBD<-convertIBD(c(0.5,0.5,0.0))
}
