#' Get mapping between reference and gene using TranscriptDb  
#' 
#' Gets the mapping between reference (usually transcripts) and gene's EntrezID using a TranscriptDb object.
#'
#' @param txDb A \code{TranscriptDb} object.
#' @param exons A logical value, indicating if exons should be extracted.
#' @param verbose A logical value, indicating if messages should be printed. Default is \code{FALSE}.
#' @return
#' \itemize{
#' \item df_GENETR A \code{data.frame} containing the genomic ranges of the transcripts, ordered by chromosome and strand;
#' \item df_EXON A \code{data.frame} containing the genomic ranges of the exons, ordered by chromosome and strand.
#' \item stack_EXON A \code{data.frame} containing mappings between exons in \code{df_EXON} and transcripts in \code{df_GENETR}.
#' \item txDb_Key A character vector, giving the EntrezID's column name.
#' \item txGroup A \code{data.frame} containing the index of each chromosome and strand pair \code{df_GENETR}. Used to speed up computation 
#' }
#' @author Alexandre Thiery
#' @keywords internal
#' @import GenomicFeatures
#' @import AnnotationDbi
#' @import RSQLite
annotateGetRefFromTranscriptDb <- function(txDb, exons, verbose = FALSE)
{
    if(verbose) message("Get reference data...")
    if(missing(txDb))
    {
        suppressPackageStartupMessages(expr = library(TxDb.Hsapiens.UCSC.hg19.knownGene))
        txDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    if(!inherits(x = txDb, what = "TranscriptDb"))
    {
        stop("'txDb' is not in the correct format.'")
    }
    ## get GRanges
    transcriptsTemp <- NULL
    txDb_Key = "GENEID"
    COLS <- c('GENEID', 'TXNAME')
    if(exons) COLS <- c(COLS, 'EXONID')
    transcriptsTemp <- GenomicFeatures::transcripts(x = txDb, columns = COLS)
    ## remove exotic chromosomes 
    temp = grep(pattern = "_", x = as.character(transcriptsTemp@seqnames), fixed = TRUE)
    if(length(temp) > 0) {
        transcriptsTemp <- transcriptsTemp[-temp,]
    }
    
    ## to data.frame 'transcriptsTemp'
    df_GENETR <- GenomicRanges::as.data.frame(transcriptsTemp)
    colnames(df_GENETR)[colnames(df_GENETR) == "seqnames"] = "chr"
    ## set correct format
    if(!is.character(df_GENETR$chr)) df_GENETR$chr <- as.character(df_GENETR$chr)
    if(!is.character(df_GENETR$strand)) df_GENETR$strand <- as.character(df_GENETR$strand)
    if(!is.character(df_GENETR$GENEID)) df_GENETR$GENEID <- as.character(df_GENETR$GENEID) 
    if(!is.character(df_GENETR$TXNAME)) df_GENETR$TXNAME <- as.character(df_GENETR$TXNAME)
    if(!is.integer(df_GENETR$start)) df_GENETR$start <- as.integer(df_GENETR$start)
    if(!is.integer(df_GENETR$end)) df_GENETR$end <- as.integer(df_GENETR$end) 
    
    ## order by chromosome and strand
    df_GENETR$txGroup <- paste(df_GENETR$chr, df_GENETR$strand, sep = "")
    if(!is.character(df_GENETR$txGroup)) df_GENETR$txGroup <- as.character(df_GENETR$txGroup)
    df_GENETR <- df_GENETR[order(df_GENETR$txGroup),]
    
    if(exons)
    {
        exonsTemp <- GenomicFeatures::exons(x = txDb, columns = c('EXONID'))
       ## to data.frame 'exonsTemp'
        temp = grep(pattern = "_", x = as.character(exonsTemp@seqnames), fixed = TRUE)
        if(length(temp) > 0) {
            exonsTemp <- exonsTemp[-temp,]
        }
        df_EXON = as.data.frame(exonsTemp)
        rm(exonsTemp)
        df_EXON = df_EXON[,c("seqnames","start","end","EXONID")]
        colnames(df_EXON)[1] = "chr"
        ## to data.frame 'exonsTemp'
       
       stack_EXON = GenomicFeatures::as.list(transcriptsTemp$EXONID)
       ##stack_EXON[[2]]
       ##stack_EXON = as.data.frame(unstack(transcriptsTemp$EXONID))
       ##stack_EXON$name = df_GENETR$TXNAME[ as.numeric( stack_EXON$name) ]
        
        ##if(!is.character(stack_EXON$name)) stack_EXON$name <- as.character(stack_EXON$name)
        ##if(!is.integer(stack_EXON$value)) stack_EXON$value <- as.integer(stack_EXON$value)
        ##
        if(!is.character(df_EXON$EXONID)) df_EXON$EXONID <- as.character(df_EXON$EXONID)
        if(!is.integer(df_EXON$start)) df_EXON$start <- df_EXON(df_GENETR$start)
        if(!is.integer(df_EXON$end)) df_EXON$end <- df_EXON(df_GENETR$end) 
    } else {
        df_EXON <- data.frame(chr = c(), start = c(), end = c(), EXONID = c())
        ##stack_EXON <- data.frame(value = c(), name = c())
        stack_EXON <- vector(mode = "list", length = 0)
    }
    rm(transcriptsTemp)
    
   
    txGroupIndex = sapply(unique(df_GENETR$txGroup), function(x) range(which(df_GENETR$txGroup == x)))-1
    txGroup = data.frame(group = colnames(txGroupIndex), index1 = txGroupIndex[1,], index2 = txGroupIndex[2,])
    
    return(list(df_GENETR = df_GENETR, txDb_Key = txDb_Key, txGroup = txGroup, df_EXON = df_EXON, stack_EXON = stack_EXON))
}