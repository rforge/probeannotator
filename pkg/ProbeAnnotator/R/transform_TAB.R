#' @title transform_TAB
#'
#' A function modify platforms' TAB files to FASTA
#'
#' @name transform_TAB
#' @param fasta.file a character vector giving the name of the input TAB file.
#' @param file.out a character vector giving the name of the ouput FASTA file.
#' @param probe_column a character vector giving the names of the probe column to extract. 
#' @param seq_column a character vector giving the names of the sequence column to extract. 
#' @param hasProbeset a logical value indicating if the platform is arranged in probesets or not. 
#' @details
#' The aim of function step is to transform the platforms sequence TAB files in fasta files format. For example, taking a sequence name from Mouse:
#' @export
#' @useDynLib ProbeAnnotator
transform_TAB = function(fasta.file, file.out, probe_column, seq_column, hasProbeset) {
  # Get the number of transcripts ...
  message("Transform TAB file")

  # Get header
  message(" - get header information")
  header = readLines(fasta.file,n = 1)
  columns = strsplit(header,split="\t",fixed=TRUE)[[1]]
  
  
  #header = paste("ALEX","B1","B2","B3",sep = "\t")
  #probe_column = "B1"
  #seq_column = "B3"
  
  # Get pattern from desier probe and seq column
  probe_index = which(columns == probe_column)
  seq_index = which(columns == seq_column)
  pattern = rep("%*s",length=length(columns))
  pattern[c(probe_index,seq_index)] <- "%s"
  pattern <- paste(pattern,collapse="\t")
  
  
  probe_first = 0
  if(probe_index > seq_index)
    probe_first=1
  
  probe_set = 0 #no probe_set
  if(hasProbeset)
    probe_set = 1 #has probe_set
  
  args <- c(probe_first, probe_set)
  
  message(" - print ouput file")
  res3 <- .C("transform_TAB", as.character(fasta.file), as.character(file.out), 
             as.character(pattern),as.integer(args))
}