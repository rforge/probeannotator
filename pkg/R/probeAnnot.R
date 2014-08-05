#' @title transform_UCSC
#'
#' A function to annotate probes from different platforms using results from Bowtie.
#'
#' @name probeAnnot
#' @param fasta.file a character vector giving the name of the probes' FASTA file.
#' @param bowtie.files a character vector giving the names of the probes' BOWTIE alignement map files.
#' @param file.out a character vector giving the name of the ouput file.
#' @param probe_pattern a character vector giving the PCRE expression for the probes (used only if you have probesets).
#' @param gene_pattern a character vector giving the PCRE expression for the genes.
#' @param weights a vector of double giving the weights associated to the number of missmatch (values must be between 0 and 1).
#' @export
#' @useDynLib ProbeAnnotator
probeAnnot = function(fasta.file, bowtie.files, file.out, probe_pattern = "", gene_pattern = "", weights = c(1,0.5,0.25)) {
  #PARAMETERS
  ## - FASTA FILE
  fasta.file = path.expand(fasta.file)
  if(!file.exists(fasta.file)) {
    message(sprintf("Cannot find fasta file '%s'.\n", fasta.file))
    stop()
  }
  ## - BOWTIE FILE
  bowtie.files = path.expand(bowtie.files)
  if(sum(!file.exists(bowtie.files)) > 0) {
    message(sprintf("Cannot find bowtie file(s) '%s'.\n", bowtie.files[!file.exists(bowtie.files)]))
    stop()
  }
  ## - FILE OUT
  file.out = path.expand(file.out)
  
  ## - PROBE PATTERN
  if(probe_pattern == "") {
    # if 'probe_pattern' is empty string, then
    # it means platform has no probeset 
    # (i.e. each probe will be treated individually)
    hasP = 0
  } else {
    # otherwise, the platform has probesets
    hasP = 1
  }
  ## - GENE PATTERN
  if(gene_pattern == "")
  {
    # if 'gene_pattern' is empty string, then set
    # pattern to match the entire character string
    gene_pattern = "^.+$"
  }
  ## - WEIGHTS
  # stop if:
  # - some weights are negative
  # - there are zero weigths
  # - all are equal to 0
  # change if:
  # - the max weigths is higher than 1 (normalise to 1)
  rw = range(weights)
  maxMissmatch = length(weights)
  if(maxMissmatch == 0) { 
    message("The weights vector is of length 0.\n")
    stop()
  }
  if(sum(weights < 0) > 0) {
    message("Error: Some weights are negative.\n")
    stop()
  }
  if(sum(weights) == 0) {
    message("The weights are all equal to 0.\n")
    stop()
  }
  if(max(rw) > 1)
  {
    warning("Weights' maximum is greater than 1. Normalising weights.")
    # rw[rw < 0] = 0
    # rw = rw/max(rw) 
    weights <- weights/max(rw)
  }
  ## - Number of BOWTIE files 'nb_map'
  nb_map = length(bowtie.files)
  
  #FUNCTION CALL
  res <- .C("score", as.character(fasta.file), as.character(bowtie.files), as.integer(nb_map),
            as.character(gene_pattern), as.character(probe_pattern), as.integer(hasP), 
            as.double(weights), as.integer(maxMissmatch),
            as.character(file.out))    
}