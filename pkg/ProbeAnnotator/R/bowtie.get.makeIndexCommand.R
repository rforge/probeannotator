#' @title bowtie.get.makeIndexCommand
#'
#' A function to create the command to create an index with Bowtie.
#'
#' @name bowtie.get.makeIndexCommand
#' @param file.in a character vector giving the name of the input FASTA file.
#' @param index a character vector giving the name of BOWTIE index to create.
#' @param bowtie.path a character vector giving the path of BOWTIE source files. 
#' @param index.dir a character vector giving the path of BOWTIE indexes files. Default value is \code{bowtie.path + "/indexes/"}
bowtie.get.makeIndexCommand <- function(file.in, index, bowtie.path, index.dir = NULL) {
  
  #parameters
  if(is.null(index.dir))
    index.dir <- file.path(bowtie.path, "indexes")

  #checks
  index.dir <- tryNormalizePath(index.dir,FALSE)
  bowtie.path <- tryNormalizePath(bowtie.path,FALSE)
  file.in <- tryNormalizePath(file.in,TRUE)
  
  error <- any(sapply(list(index.dir, bowtie.path, file.in), is.null))
  if(error) stop()
  
  #build command
  cmd.index   = file.path(index.dir, index)
  cmd.bowtie  = file.path(bowtie.path,"bowtie-build")
  cmd = paste(cmd.bowtie, "-f", file.in , cmd.index, sep = " ")
  return(cmd)
}
