\name{make.profile}
\alias{make.profile}
\alias{make.profiles}
\title{Convert genotype profile to reference profile format}
\description{Convert genotype profile(s) to reference profile format}
\usage{
make.profile(gt,name='K')
make.profiles(typed.gts)
}
\arguments{
  \item{gt}{genotype profile}
  \item{name}{character string used to name profile in output data frame}
  \item{typed.gts}{named list of genotype profiles}
}
\value{data frame containing reference profile(s)}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\examples{
data(test2data)
S1prof<-make.profile(S1gt,'S1')
C<-0.001
mixD<-DNAmixture(list(epg),k=3,K='S1',reference.profile=S1prof,C=list(C),database=db)
}
