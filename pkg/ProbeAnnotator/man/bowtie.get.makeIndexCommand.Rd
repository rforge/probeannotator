\name{bowtie.get.makeIndexCommand}
\alias{bowtie.get.makeIndexCommand}
\title{bowtie.get.makeIndexCommand

A function to create the command to create an index with Bowtie.}
\usage{
bowtie.get.makeIndexCommand(file.in, index, bowtie.path, index.dir = NULL)
}
\arguments{
  \item{file.in}{a character vector giving the name of the
  input FASTA file.}

  \item{index}{a character vector giving the name of BOWTIE
  index to create.}

  \item{bowtie.path}{a character vector giving the path of
  BOWTIE source files.}

  \item{index.dir}{a character vector giving the path of
  BOWTIE indexes files. Default value is \code{bowtie.path
  + "/indexes/"}}
}
\description{
bowtie.get.makeIndexCommand

A function to create the command to create an index with
Bowtie.
}

