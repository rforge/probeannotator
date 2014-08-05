#' @title bowtie.get.makeAlignmentCommand
#'
#' A function to create the command to perform alignement with Bowtie.
#'
#' @name bowtie.get.makeAlignmentCommand
#' @param file.in a character vector giving the name of the input FASTA file.
#' @param file.out a character vector giving the name of the output BOWTIE map file.
#' @param index a character vector giving the name of BOWTIE index to align with.
#' @param bowtie.path a character vector giving the path of BOWTIE source files.
#' @param index.dir a character vector giving the path of BOWTIE indexes files. Default value is \code{bowtie.path + "/indexes/"}
#' @param p a numerical value giving the number of core used for multithreading. If \code{p=0}, then multithreading is not used. Default is code{p=0}.
bowtie.get.makeAlignmentCommand <- function(file.in, file.out, index, bowtie.path, index.dir = NULL,p=0) {
  #parameters
  if(is.null(index.dir))
    index.dir <- file.path(bowtie.path, "indexes")

  #checks 
  index.dir <- tryNormalizePath(index.dir,FALSE)
  bowtie.path <- tryNormalizePath(bowtie.path,FALSE)
  file.in <- tryNormalizePath(file.in,TRUE)
  file.out <- tryNormalizePath(file.out,TRUE,FALSE)
  
  error <- any(sapply(list(index.dir, bowtie.path, file.in, file.out), is.null))
  if(error) stop()
  
  #build command
  cmd.options = "-v 3 -a --best"
  if(p > 0)
    cmd.options = paste(cmd.options, " -p ", p)
  cmd.bowtie  = file.path(bowtie.path,"bowtie")
  cmd.index   = file.path(index.dir, index)
  cmd     = paste(cmd.bowtie, cmd.options, cmd.index, "-f", file.in, file.out, sep = " ")
  return(cmd)
}
