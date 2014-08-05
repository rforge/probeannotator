\name{transform_UCSC}
\alias{transform_UCSC}
\title{transform_TAB

A function modify platforms' TAB files to FASTA}
\usage{
transform_TAB(fasta.file, file.out, probe_column, seq_column, hasProbeset)

transform_UCSC(fasta.file, file.out, appendString = "", organism = "Mm")
}
\arguments{
  \item{fasta.file}{a character vector giving the name of
  the input TAB file.}

  \item{file.out}{a character vector giving the name of the
  ouput FASTA file.}

  \item{probe_column}{a character vector giving the names
  of the probe column to extract.}

  \item{seq_column}{a character vector giving the names of
  the sequence column to extract.}

  \item{hasProbeset}{a logical value indicating if the
  platform is arranged in probesets or not.}

  \item{fasta.file}{a character vector giving the name of
  the input FASTA file.}

  \item{file.out}{a character vector giving the name of the
  ouput FASTA file.}

  \item{appendString}{a character vector giving the string
  to append to each refSeq and gene_id.}

  \item{organism}{a character vector giving the name of the
  organism: human (Hs) ou mouse (Mm). Default is
  \code{Mn}.}
}
\description{
transform_TAB

A function modify platforms' TAB files to FASTA

transform_UCSC

A function modify Fasta files obtained from UCSC
}
\details{
The aim of function step is to transform the platforms
sequence TAB files in fasta files format. For example,
taking a sequence name from Mouse:

The aim of function step is to transform the sequences
names in fasta files downloaded from UCSC. For example,
taking a sequence name from Mouse: convert:
mm10_refGene_NM_001282945; to : NM_001282945|11539. That,
the sequences' names will be \code{refGene|GeneID}.
}

