\name{bowtie.makeAlignment}
\alias{bowtie.makeAlignment}
\title{bowtie.makeAlignment

A function to create and execute bowtie alignment.}
\usage{
bowtie.makeAlignment(file.in, file.out, index,
  bowtie.path = "/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0",
  index.dir = NULL, p = 0)
}
\arguments{
  \item{file.in}{a character vector giving the name of the
  input FASTA file.}

  \item{file.out}{a character vector giving the name of the
  output BOWTIE map file.}

  \item{index}{a character vector giving the name of BOWTIE
  index to align with.}

  \item{bowtie.path}{a character vector giving the path of
  BOWTIE source files. Default value for the BCF is
  \code{"/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0"}}

  \item{index.dir}{a character vector giving the path of
  BOWTIE indexes files. Default value for the BCF is
  \code{bowtie.path + "/indexes/"}.}

  \item{p}{a numerical value giving the number of core used
  for multithreading. If \code{p=0}, then multithreading is
  not used. Default is code{p=0}.}
}
\description{
bowtie.makeAlignment

A function to create and execute bowtie alignment.
}

