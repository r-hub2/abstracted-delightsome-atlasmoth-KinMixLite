\name{checkdata}
\alias{checkdata}
\title{Check data for absences of markers or allele values}
\description{Check input data used by KinMix for absences of required markers or allele values.}
\usage{
checkdata(epg,database,typed.gts=NULL,reference.profiles=NULL)
}
\arguments{
  \item{epg}{data frame, the epg; see \code{\link{formats}}.}
  \item{database}{data frame, the db; see \code{\link{formats}}.}
  \item{typed.gts}{named list of genotype objects; see \code{\link{formats}.}}
  \item{reference.profiles}{data frame containing reference.profiles in DNAmixturesLite format}
}
\value{NULL}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\examples{
data(test2data)

checkdata(epg,db,list(C=Cgt,F=Fgt),make.profiles(list(M=Mgt)))
}
