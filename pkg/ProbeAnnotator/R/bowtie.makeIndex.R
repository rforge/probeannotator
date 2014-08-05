#' @title bowtie.makeIndex
#'
#' A function to create and execute bowtie index.
#'
#' @name bowtie.makeIndex
#' @param file.in a character vector giving the name of the input FASTA file.
#' @param index a character vector giving the name of BOWTIE index to align with.
#' @param bowtie.path a character vector giving the path of BOWTIE source files.
#' @param index.dir a character vector giving the path of BOWTIE indexes files. Default value is \code{bowtie.path + "/indexes/"}
#' @export
bowtie.makeIndex <- function(file.in, index, bowtie.path, index.dir = NULL) {
  cmd = bowtie.get.makeIndexCommand(file.in, index, bowtie.path, index.dir)
  if(is.null(cmd))
    stop()
  
  
  system(cmd)
}
