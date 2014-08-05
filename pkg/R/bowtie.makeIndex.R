#' @title bowtie.makeIndex
#'
#' A function to create and execute bowtie index.
#'
#' @name bowtie.makeIndex
#' @param file.in a character vector giving the name of the input FASTA file.
#' @param index a character vector giving the name of BOWTIE index to align with.
#' @param bowtie.path a character vector giving the path of BOWTIE source files. Default value for the BCF is \code{"/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0"}
#' @param index.dir a character vector giving the path of BOWTIE indexes files. Default value for the BCF is \code{bowtie.path + "/indexes/"}
#' @export
bowtie.makeIndex <- function(file.in, index, bowtie.path = "/net/sib-pc12/export/big/csoneson/Projects/almac_annotation/bowtie-1.0.0", index.dir = NULL) {
  cmd = bowtie.get.makeIndexCommand(file.in, index, bowtie.path, index.dir)
  if(is.null(cmd))
    stop()
  
  
  system(cmd)
}
