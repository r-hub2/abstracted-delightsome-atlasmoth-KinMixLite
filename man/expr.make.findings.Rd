\name{expr.make.findings}
\alias{expr.make.findings}
\title{Coding additional findings as expression}
\description{Returns an expression that will be evaluated 
in logL.UKX whenever the likelihood of the model is calculated using the
current method, and encodes the additional findings needed to implement the 
method; the details of the model and the extra information needed are held in the list z
}
\usage{
expr.make.findings(z)
}
\arguments{
  \item{z}{A list specifying the additional findings; for the format, see Details}
}
\details{
Each component of the list \code{z} is a list encoding a particular type of additional finding: the first component of this (sub-)list being a character string specifying the type of finding, and the remainder of its components being named parameters giving details of the finding. The types of finding and the valid parameters of each are as follows:
\describe{
\item{\code{Male}}{\code{ind}: index of relevant contributor: which `unknown' contributor are we modelling by amending his CPTs?}
\item{\code{Female}}{\code{ind}: index of relevant contributor}
\item{\code{Rlikd}}{\code{aca}: allele count array, \code{cgt}: character string naming genotype profile data frame, \code{evid}: character string naming list with one component for each marker, whose value is the evidence}
\item{\code{Aca}}{\code{ind}: index of relevant contributor, \code{aca}: allele count array}
\item{\code{Caca}}{\code{ind}: index of relevant contributor, \code{aca}: allele count array}
\item{\code{Denom}}{no parameters}
}
If \code{z} is NULL, then there are no additional findings.
}
\value{Expression encoding the additional findings.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
