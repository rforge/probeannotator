#' Fast Gene Annotation by Chromosomal Location
#' @description This method allows fast gene annotation from chromosomal location
#' 
#' @param x A \code{GRanges} object or \code{data.frame}, giving the coordinates to annotate.
#' @param txDb A \code{TranscriptDb} object, giving the genomic references for the alignement. If \code{txDb} is missing, the default \code{TxDb.Hsapiens.UCSC.hg19.knownGene} will be used.
#' @param orgDb A \code{OrganismDb} object, giving the details of the genomic references. If \code{orgDb} is missing, the default \code{org.Hs.eg.db} will be used.
#' @param orgDb_Columns A character vector (optional), giving which columns to extract from the \code{orgDb}. Note that if \code{orgDb_Columns} is used, then the user selected columns will be selected instead of the default request on \code{org.Hs.eg.db}, see details section.
#' @param txAnnotRange A integer vector, giving the window size for the genes' promotor and extended sites. Default is \code{txAnnotRange= c(1500,2000)}, See details.
#' @export
#' @import GenomicFeatures
#' @importFrom GenomicRanges as.data.frame
#' @import Rcpp
#' @useDynLib ProbeAnnotator
annotateByLocation = function(x, txDb, txAnnotRange = c(1500,2000), orgDb, orgDb_Columns)
{
    sep_intra = ";"
    sep_inter = '\\\\'
    
    ## checking x
    if(missing(x)) stop("'x' is missing.")
    if(inherits(x = x, what = "GRanges"))
    {
        x <- GenomicRanges::as.data.frame(x)
        colnames(x)[colnames(x) == "seqnames"] = "chr"
        if("ID" %in% colnames(x)) stop("No column 'ID' in 'x'.")
    } else if(is.data.frame(x))
    {
        temp = c("ID", 'chr', 'strand', 'probeEnd', 'start' ,'end') %in% colnames(x)
        if(any(!temp)) 
        {
            stop(sprintf("Column(s) %s not found in 'x'", paste( colnames(x)[!temp] , collapse = ", " )))
        }
    }
    else {
        stop("'x' is not in the correct format.")
    }
    x$txGroup = paste(x$chr, x$strand, sep = "")
    
    message("Get reference data...")
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
    txGRanges <- NULL
    txDb_Key = "GENEID"
    txGRanges <- genes(x = txDb, columns = c('GENEID'))
   
    ## to data.frame
    txGRanges <- GenomicRanges::as.data.frame(txGRanges)
    colnames(txGRanges)[colnames(txGRanges) == "seqnames"] = "chr"
    ## remove exotic chromosomes 
    temp = grep(pattern = "_", x = txGRanges$chr, fixed = TRUE)
    if(length(temp) > 0) {
        txGRanges <- txGRanges[-temp,]
    }
    ## order by chromosome and strand
    txGRanges$txGroup <- paste(txGRanges$chr, txGRanges$strand, sep = "")
    txGRanges <- txGRanges[order(txGRanges$txGroup, txGRanges$start),]
    txGroupIndex = sapply(unique(txGRanges$txGroup), function(x) range(which(txGRanges$txGroup == x)))-1
    txGroup = data.frame(group = colnames(txGroupIndex), index1 = txGroupIndex[1,], index2 = txGroupIndex[2,])
    ## set correct format
    if(!is.character(txGRanges$chr)) txGRanges$chr <- as.character(txGRanges$chr)
    if(!is.character(txGRanges$strand)) txGRanges$strand <- as.character(txGRanges$strand)
    if(!is.character(txGRanges[,txDb_Key])) txGRanges[,txDb_Key] <- as.character(txGRanges[,txDb_Key]) 
    if(!is.integer(txGRanges$start)) txGRanges$start <- as.integer(txGRanges$start)
    if(!is.integer(txGRanges$end)) txGRanges$end <- as.integer(txGRanges$end) 
    
    orgDb_Key = ""
    orgDb_Columns = c("")
    ## 
    message("Get organism data...")
    if(missing(orgDb))
    {
        # default query using 'org.Hs.eg.db'
        # faster than using select
        suppressPackageStartupMessages(expr = library(org.Hs.eg.db))
        orgDb <- org.Hs.eg.db
        txInfo = dbGetQuery(conn = org.Hs.eg.db::org.Hs.eg_dbconn(),
                            statement = sprintf("SELECT genes.gene_id AS 'ENTREZID',gene_name AS 'GENENAME',symbol AS 'SYMBOL',group_concat(alias_symbol,'%s') AS 'ALIAS' FROM genes JOIN gene_info ON genes._id = gene_info._id JOIN alias ON genes._id = alias._id GROUP BY genes._id", sep_intra))
        orgDb_Key = "ENTREZID"
        orgDb_Columns = c("ENTREZID", "GENENAME", "SYMBOL", "ALIAS")
        txInfo[1,]
        for(i in 1:ncol(txInfo)) txInfo[,i] <- as.character(txInfo[,i])
    } else
    {
        if(!inherits(x = orgDb, what = "OrgDb"))
        {
            stop("'OrgDb' is not in the correct format.")
        } 
        if(missing(orgDb_Columns)) 
        {
            stop("'orgDb_Columns' is missing.")
        } 
        if (is.character(orgDb_Columns))
        {
            ## check if columns are in orgDb
            temp  = orgDb_Columns %in% columns(x = orgDb)
            if(any(!temp)) 
            {
                stop(sprintf("Column(s) %s not found in 'orgDb'", paste( orgDb_Columns[!temp] , collapse = ", " )))
            }
        } else {
            stop("'orgDb_Columns' is not in the correct format.")
        }
        
        orgDb_Key = "ENTREZID"
        if(!(orgDb_Key %in% orgDb_Columns))
        {
            orgDb_Columns <- c(orgDb_Key, orgDb_Columns)
        } else { 
            #set orgDb_Key in first pos
            temp1 = which(orgDb_Columns == orgDb_Key)
            temp2 = orgDb_Columns[1]
            orgDb_Columns[1] = orgDb_Key
            orgDb_Columns[temp1] = temp2
        }
        txInfo = select(orgDb, keys(orgDb, orgDb_Key), 
                        columns = orgDb_Columns, 
                        keytype = orgDb_Key)
        txInfo = as.data.frame(
            do.call("rbind", 
                    lapply( unique(txInfo[,1]), function(y) {
                        which1 = which(txInfo[,1] == y)
                        c(y, apply(txInfo[which1,-1], 2, function(z) paste(unique(z), collapse = sep_intra)))
                    })))
        colnames(txInfo)[1] = orgDb_Key
    }
    
    message("Perform annotation...")
    result <- .Call('ProbeAnnotator_annotateByLocation', PACKAGE = 'ProbeAnnotator',
                    x, txGRanges, txGroup, txDb_Key, txAnnotRange, 
                        txInfo, orgDb_Key, orgDb_Columns)
    message("Done")
    as.data.frame(do.call("cbind", result))    
}