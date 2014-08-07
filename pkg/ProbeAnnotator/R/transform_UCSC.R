#' @title transform_UCSC
#'
#' A function modify Fasta files obtained from UCSC
#'
#' @name transform_UCSC
#' @param fasta.file a character vector giving the name of the input FASTA file.
#' @param file.out a character vector giving the name of the ouput FASTA file.
#' @param appendString a character vector giving the string to append to each refSeq and gene_id.
#' @param organism a character vector giving the name of the organism: human (Hs) ou mouse (Mm). Default is \code{Mn}.
#' @details
#' The aim of function step is to transform the sequences names in fasta files downloaded from UCSC. For example, taking a sequence name from Mouse:
#' convert: mm10_refGene_NM_001282945;
#' to     : NM_001282945|11539.
#' That, the sequences' names will be \code{refGene|GeneID}.
#' @export
#' @importFrom AnnotationDbi as.list
#' @useDynLib ProbeAnnotator
transform_UCSC = function(fasta.file, file.out, appendString = "", organism = "Mm") {


  if(!TRUE)
	stop("PCRE lib is not available")

  # Get the number of transcripts ...
  message("Extract the refSeqs count")
  result_length = 0
  
  if(.Platform$OS.type == "unix") {
    cmd1 = sub("FILE", fasta.file, "tr -d -c '>' < FILE | awk '{ print length; }'")
    result_length = as.integer(system(cmd1,intern=TRUE ))
  } else {
    res1 = .C("get_UCSC_refSeq_count", as.character(fasta.file), as.integer(result_length))
    result_length = res1[[2]] 
  }
  
  empty_res = paste(collapse = "", rep("x",30))
  result = rep(empty_res  , result_length)
  
  # Get the transcripts ...
  message("Extract the refSeqs' name")
  res2 <- .C("get_UCSC_refSeq", as.character(fasta.file),
             as.character(result), as.integer(result_length))
  
  
  # Get their GENE_ID
  x <- NULL
  message("Get refSeqs' entrez id")
  if(organism == "Mm") {
    require(org.Mm.eg.db,warn.conflicts=FALSE, quietly=TRUE)
    x <- org.Mm.egREFSEQ2EG
  } else if(organism == "Hs") {
    require(org.Hs.eg.db,warn.conflicts=FALSE, quietly=TRUE)
    x <- org.Hs.egREFSEQ2EG
  }
  #sql2 = "SELECT gene_id, accession FROM genes JOIN refseq ON refseq._id = genes._id"
  #sel2 = dbGetQuery(org.Hs.eg.db::org.Hs.eg_dbconn(), sql2)

  # Get the RefSeq identifier that are mapped to an entrez gene ID
  mapped_seqs <- mappedkeys(x)
  # Convert to a list
  # AnnotationDbi:::as.list
  xx <- as.list(x[mapped_seqs])
  GENE_ID = xx[res2[[2]]]
  REFS_ID = res2[[2]]
  
  #rm(xx, res2, result);
  wn = which(is.na(names(GENE_ID)))
  GENE_ID[wn] = paste("unknown_", 1:length(wn), sep = "")
  
  u_rep = unique(REFS_ID[duplicated(REFS_ID)])
  for(i in u_rep) {
    w = which(REFS_ID == i)
    REFS_ID[w]  = paste(i, ".", 1:length(w), sep = "")
  }
  
  
  final_result = c()
  # Get their new id : refSeq|GENE_ID
  if(appendString == "") {
    final_result = paste(REFS_ID,GENE_ID,sep = "|")
  } else {
    final_result = paste(REFS_ID, "_", appendString, "|", GENE_ID, "_", appendString,sep = "" )
  }
  grep1 = grep("(", final_result, fixed = TRUE)
  if(length(grep1) > 0 ) {
    f.temp <- sapply(grep1, function(x) {
      paste(REFS_ID[x], "_", appendString, "|", GENE_ID[[x]][1], "_", appendString,sep = "" )
    })
    final_result[grep1] = f.temp
  }
 
  # print result
  message("Print ouput file")
  res3 <- .C("print_transformed_UCSC", as.character(fasta.file), as.character(file.out), 
             as.character(final_result))
}
