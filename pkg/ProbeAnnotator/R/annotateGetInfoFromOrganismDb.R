#' Get gene information using EntrezID and OrganismDb  
#' 
#' Gets the genes information in an OrganismDb using EntrezID keys.
#'
#' @param orgDb A \code{OrganismDb} object. If \code{orgDb} is missing, the default \code{org.Hs.eg.db} will be used.
#' @param orgDb_Columns A character vector, giving which columns to extract from the \code{orgDb}. Note that if \code{orgDb_Columns} is used, then the user selected columns will be selected instead of the default request on \code{org.Hs.eg.db}.
#' @param sep_intra A character vector, giving the separator character for gene information, see details.
#' @param verbose A logical value, indicating if messages should be printed. Default is \code{FALSE}.
#' @return
#' \itemize{
#' \item txInfo A \code{data.frame} containing the genes' information.
#' \item orgDb_Key A character vector, giving the EntrezID's column name.
#' \item orgDb_Columns A character vector, giving the selected column names (see Arguments section).
#' }
#' @author Alexandre Thiery
#' @keywords internal
#' @import GenomicRanges
annotateGetInfoFromOrganismtDb <- function(orgDb, orgDb_Columns, sep_intra, verbose = FALSE) {
    orgDb_Key = ""
    ## 
    if(verbose) message("Get organism data...")
    if(missing(orgDb))
    {
        # default query using 'org.Hs.eg.db'
        # faster than using select
        suppressPackageStartupMessages(expr = library(org.Hs.eg.db))
        orgDb <- org.Hs.eg.db
        txInfo = RSQLite::dbGetQuery(conn = org.Hs.eg.db::org.Hs.eg_dbconn(),
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
        txInfo = AnnotationDbi::select(orgDb, keys(orgDb, orgDb_Key), 
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
    return(list(txInfo = txInfo, orgDb_Key = orgDb_Key, orgDb_Columns = orgDb_Columns))
}