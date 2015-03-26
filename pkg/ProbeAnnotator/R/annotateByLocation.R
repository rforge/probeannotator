#' @rdname annotateBy
#' @examples
#' ## Example of 3 coordinates of the MGMT gene on chromosome 10
#' ## - coordinates: chr10:27132612-2713562
#' ## - assembly: hg19 (default in annotateByLocation)
#' start = c(131263500, 131264960, 131265460)
#' probeID = sprintf("probe_%d",1:3)
#' 
#' ## Using x in GRanges format
#' gr = GRanges(seqnames = "chr10",  
#'              strand = "+",
#'              ranges = IRanges(start = start,
#'                               width = 20),
#'              ID = probeID)
#' annot_gr = annotateByLocation(x = gr, mapType = "EXON")
#' 
## Using x in data.frame format
#' df = data.frame(chr = "chr10", 
#'                 strand = "+", 
#'                 start = start, 
#'                 end = start+20, 
#'                 ID = probeID)
#' annot_df = annotateByLocation(x = df, mapType = "EXON")
#' 
#' ## Check if both results are the same
#' all.equal(annot_gr, annot_df)
#' 
#' print(annot_gr)
#' @export
#' @import GenomicFeatures
#' @importFrom GenomicRanges as.data.frame
#' @import Rcpp
#' @useDynLib ProbeAnnotator
annotateByLocation = function(x, txDb, mapType, promotorRange = 1500, extendedRange = 2000, orgDb, orgDb_Columns, sep_intra = ";", sep_inter = '\\', verbose = FALSE)
{
    ## check map
    mapType = .checkVector(mapType, "is.character", name = quote(mapType))
    maps = c("EXONINTRON",  "NO_EXONINTRON", "EXON")
    mapType = match.arg(arg = mapType, choices = maps, several.ok = FALSE)
    mapTypeIndex = which(mapType == maps)-1
    exons = mapType %in% c("EXONINTRON", "EXON")
    ## check promotor & extended
    promotorRange = .checkVector(promotorRange, "is.integer", name = quote(promotorRange))
    extendedRange = .checkVector(extendedRange, "is.integer", name = quote(extendedRange))
    ## check sep_intra, sep_inter
    sep_intra = .checkVector(sep_intra, "is.character", name = quote(sep_intra))
    sep_inter = .checkVector(sep_inter, "is.character", name = quote(sep_inter))
    ## check verbose
    verbose = .checkVector(verbose, "is.logical", name = quote(verbose))
    
    ## checking x
    if(missing(x)) stop("'x' is missing.")
    if(inherits(x = x, what = "GRanges"))
    {
        x <- GenomicRanges::as.data.frame(x)
        colnames(x)[colnames(x) == "seqnames"] = "chr"
        if(!("ID" %in% colnames(x))) stop("No column 'ID' in 'x'.")
    } else if(is.data.frame(x))
    {
        temp = c("ID", 'chr', 'strand', 'end', 'start') %in% colnames(x)
        if(any(!temp)) 
        {
            stop(sprintf("Column(s) %s not found in 'x'", paste( c("ID", 'chr', 'strand', 'end', 'start')[!temp] , collapse = ", " )))
        }
    }
    else {
        stop("'x' is not in the correct format.")
    }
    ## set correct format
    if(!is.character(x$chr)) x$chr <- as.character(x$chr)
    if(!is.character(x$strand)) x$strand <- as.character(x$strand)
    if(!is.character(x$ID)) x$ID <- as.character(x$ID) 
    if(!is.integer(x$start)) x$start <- as.integer(as.character(x$start)) 
    if(!is.integer(x$end)) x$end <- as.integer(as.character(x$end)) 
    x$txGroup = paste(x$chr, x$strand, sep = "")
    
    refList= annotateGetRefFromTranscriptDb(txDb, exons = exons, verbose = verbose)
    infList= annotateGetInfoFromOrganismtDb(orgDb, orgDb_Columns, sep_intra, verbose = verbose)
    
    minRange = max(c(promotorRange, extendedRange, 0))
    maxRange = max(c(extendedRange, 0))
    txAnnotRange = c(promotorRange, extendedRange, minRange, maxRange)
    txAnnotRange = as.integer(txAnnotRange)
    
    
    if(verbose) message("Perform annotation...")
    result <- .Call('ProbeAnnotator_annotateByLocation', 
                    x, refList$df_GENETR, refList$txGroup, refList$df_EXON, refList$stack_EXON, 
                    infList$txInfo, infList$orgDb_Columns,
                    c(sep_intra,sep_inter),                
                    txAnnotRange, 
                    list( c(mapTypeIndex, as.integer(verbose)) )
    )
    
    ##result <- .Call('ProbeAnnotator_annotateByLocation', PACKAGE = 'ProbeAnnotator',
    ##                x, refList$txGRanges, refList$txGroup, refList$txDb_Key, txAnnotRange, 
    ##                    infList$txInfo, infList$orgDb_Key, infList$orgDb_Columns,
    ##                sep_inter)
    if(verbose) message("Done")
    if(!is.null(result))
        return(as.data.frame(do.call("cbind", result)))    
}