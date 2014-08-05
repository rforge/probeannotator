\name{bowtie.makeIndex}
\alias{bowtie.makeIndex}
\title{bowtie.makeIndex

A function to create and execute bowtie index.}
\usage{
bowtie.makeIndex(file.in, index,
  bowtie.path = "/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0",
  index.dir = NULL)
}
\arguments{
  \item{file.in}{a character vector giving the name of the
  input FASTA file.}

  \item{index}{a character vector giving the name of BOWTIE
  index to align with.}

  \item{bowtie.path}{a character vector giving the path of
  BOWTIE source files. Default value for the BCF is
  \code{"/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0"}}

  \item{index.dir}{a character vector giving the path of
  BOWTIE indexes files. Default value for the BCF is
  \code{bowtie.path + "/indexes/"}}
}
\description{
bowtie.makeIndex

A function to create and execute bowtie index.
}

