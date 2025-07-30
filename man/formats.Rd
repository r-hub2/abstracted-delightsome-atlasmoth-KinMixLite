\name{formats}
\alias{formats}
\title{Formats}
\description{Formats for data objects in KinMix and KinMixLite}
\section{Formats}{
An \bold{allele frequency database} is a data frame containing variables \code{marker}, \code{allele} and \code{frequency} (character, numeric and numeric respectively).

A \bold{mixture profile} is a data frame containing variables \code{marker}, \code{allele} and \code{height} (character, numeric and numeric respectively).

A \bold{genotype profile} is a data frame containing variables \code{marker}, \code{allele1} and \code{allele2} (character, numeric and numeric respectively).

Examples of these 3 data formats are objects \code{db}, \code{epg} and \code{Cgt}, respectively, in \code{test2data}.

A \bold{allele count array} is an alternative format for a genotype as a named list of vectors, one for each marker. Each vector gives the number of each allele in the genotype, with the alleles listed in the order in which they appear in the \code{data} component of the relevant mixture object. 

An \bold{IBD pattern distribution} or \bold{IBD object} is a list with components \code{pr} (a numerical vector) and \code{patt} (an integer matrix with \code{nrow(patt)==length(pr)} and an even number of columns). The elements of \code{pr} are the probabilities of the IBD patterns in the corresponding rows of \code{patt}. Adjacent pairs of columns encode the genotypes of different individuals; equal elements in any row determine equality of the alleles; different elements denote independent draws from the gene pool. If the component \code{pr} is missing, functions \code{rpt.IBD} and \code{rpt.typed.relatives} assume the probabilities are equal.
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}


