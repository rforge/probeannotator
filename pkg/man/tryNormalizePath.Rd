\name{tryNormalizePath}
\alias{tryNormalizePath}
\title{tryNormalizePath

A function check that files or directories exists, and normalize their path.}
\usage{
tryNormalizePath(path, is.file = TRUE, check.file = TRUE)
}
\arguments{
  \item{path}{a character vector giving the name directory
  or file to check.}

  \item{is.file}{a logical value indicating wether the
  \code{path} is a directory or a file.}

  \item{check.file}{a logical value indicating (in the case
  of a file, when \code{is.file} is \code{TRUE}), if the
  file should be checked (\code{check.file=TRUE}), or if
  the directory should be checked
  (\code{check.file=FALSE})}
}
\description{
tryNormalizePath

A function check that files or directories exists, and
normalize their path.
}

