#' @title bowtie.makeAlignment
#'
#' A function to create and execute bowtie alignment.
#'
#' @name bowtie.makeAlignment
#' @param file.in a character vector giving the name of the input FASTA file.
#' @param file.out a character vector giving the name of the output BOWTIE map file.
#' @param index a character vector giving the name of BOWTIE index to align with.
#' @param bowtie.path a character vector giving the path of BOWTIE source files.
#' @param index.dir a character vector giving the path of BOWTIE indexes files. Default value is \code{bowtie.path + "/indexes/"}.
#' @param p a numerical value giving the number of core used for multithreading. If \code{p=0}, then multithreading is not used. Default is code{p=0}.
#' @export
bowtie.makeAlignment <- function(file.in, file.out, index, bowtie.path, index.dir = NULL, p=0) {
  
  cmd = bowtie.get.makeAlignmentCommand(file.in, file.out, index, bowtie.path, index.dir,p)
  if(is.null(cmd))
    stop()
  
  system(cmd)
}
