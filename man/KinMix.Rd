\name{KinMix}
\alias{KinMix}
\title{Create a DNA mixture model, with possibly related contributors}
\usage{
KinMix(data,k,C,database,K=character(0),reference.profiles = NULL, 
	contribs=NULL,typed.gts=NULL,IBD=NULL,targets=NULL,pars=NULL,mle=FALSE,
	dir=character(0),domainnamelist=NULL,
	load=FALSE,write=FALSE,dyes=NULL, 
	triangulate=TRUE,compile=TRUE,compress=TRUE,use.order=TRUE)
}
\arguments{
\item{data}{A list containing one \code{data.frame} for each DNA
mixture. Note, that in the special case of analysing just one
mixture, this still has to be specified as list(data). Each
dataset should contain variables \code{marker}, \code{allele}, and
\code{frequency}. Optionally, also a column for each reference
profile specified in \code{K}.}

\item{k}{Number of contributors.}

\item{C}{A list of thresholds, one for each mixture.}

\item{database}{A data.frame containing at least variables \code{marker}, \code{allele}, \code{frequency}.}

\item{K}{Names of reference profiles; these can be chosen freely,
but should match (possibly only a subset of) the names specified
by the reference profiles.}

\item{reference.profiles}{A data.frame containing allele counts for each reference profile, if not specified in \code{data}.}

\item{contribs}{vector of character tags identifying contributors to the mixture}

\item{typed.gts}{list of named genotype profiles}

\item{IBD}{IBD pattern distribution, or any object accepted as an argument to \code{as.IBD}}

\item{targets}{vector of character tags identifying individuals related according to \code{IBD}}

\item{pars}{optionally, a \code{mixpar} object providing parameter values for peak height model}

\item{mle}{logical: should \code{mixML} be called to estimate parameters by MLE?}

\item{dir}{Location of network files if loading or saving the networks.}

\item{domainnamelist}{Names of marker-wise network files (without hkb-extension). Default is the set of markers.}

\item{load}{Read networks from disk?}

\item{write}{Save networks as hkb files?}

\item{dyes}{A list containing a list of dyes indexed by markers}

\item{triangulate}{Triangulate the networks? Default is to
triangulate the network using a good elimination order.}

\item{compile}{Compile the networks?}

\item{compress}{Compress the network? Defaults to \code{TRUE} and is
strongly recommended for models with a large number of
contributors.}

\item{use.order}{Should the default elimination order be used for triangulation?
Otherwise the "total.weight" heuristic for triangulation in Hugin is used.}
}
\value{
An object of class \code{DNAmixture}. This contains amongst other things
\item{markers}{The joint set of markers used for the mixtures specified.}
\item{domains}{For models involving unknown contributors,
a list containing one Bayesian network (\code{hugin.domain}) per marker;
see \code{\link[DNAmixturesLite]{buildMixtureDomains}} for details on the networks}
\item{data}{A list containing for each marker the combined allele frequencies,
peak heights, and reference profiles as produced by \code{\link[DNAmixturesLite]{DNAmixtureData}}.}
It may also contain
\item{mle}{Maximum likelihood estimates of the peak height model parameters}
\item{logL}{log-likelihood of model evaluated at \code{pars}}
}
\description{
Create a DNA mixture model, with possibly related contributors.
}
\details{
Generalises \code{DNAmixture} to allow relatedness as in \code{rpt.IBD}
}
\seealso{
\code{[DNAmixturesLite]{DNAmixture}}
}
\examples{
data(test2data)

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)
baseline

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
mlD <- mixML(mixD, pars)
print(mlD$mle)
pars<-mlD$mle

baseline<-logL(mixD)(pars)
baseline

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

mixD<-KinMix(list(epg),k=2,C=list(0.001),database=db,
	contribs=c('F'),typed.gts=list(C=Cgt),IBD='parent-child',targets=c('F','C'),
	pars=pars)
log10LR<-(mixD$logL-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributor 1 is father of a typed individual Cgt 
## with mother Mgt

mixD<-KinMix(list(epg),k=2,C=list(0.001),database=db,
	contribs=c('F'),typed.gts=list(M=Mgt,C=Cgt),IBD='trio',targets=c('F','M','C'),
	pars=pars)
log10LR<-(mixD$logL-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}

