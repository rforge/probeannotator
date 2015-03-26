#' Search a "gene" contained in a character vector
#'
#' This function offers the possibility to search genes (character pattern) contained in a character vector separted by a given separator (semi-colon by default).
#' @param x A character vector.
#' @param gene A character vector, the gene (character pattern) of interest.
#' @param sep A character vector, the separator between genes (by default, \code{sep=";"}).
#' @param type A character indicating the type of values returned by the function. By default the function return a logical vector (\code{type="logical"}). The type \code{"numeric"} and \code{"value"} return position and value respectively.
#' @return The functiuon returns a vector of logical (\code{type="logical"}), numeric (\code{type="numeric"}) or character (\code{type="value"}) depending of the selected type.
#' @examples {
#' genelist <-c("MGMT","MGMT;MGMT","MGMT;EGFR","HOXA1;EGFR;HOXA2","EGFR;MGMT",NA)
#' searchGene(genelist,"EGFR",sep=";",type="logical")
#' searchGene(genelist,"EGFR",sep=";",type="numeric")
#' searchGene(genelist,"EGFR",sep=";",type="value")
#' }
#' @author Pierre Bady
#' @export
searchGene <- function(x, gene, sep=";", type = "logical") {
    type = match.arg(arg = "n", choices = c("numeric", "logical", "value"))
    gene = toupper(gene)
    x = toupper(x)
    pattern = sprintf("^%s$|^%s%s|%s%s$|%s%s%s", gene, gene, sep, sep, gene, sep, gene, sep)
    
    res = switch(type,
           numeric=grep(pattern = pattern, x = x),
           value=grep(pattern = pattern, x = x, value = TRUE),
           logical=grepl(pattern = pattern, x = x))
           
    return(res)
}
# 
# head(result$ALIAS)
# l = strsplit(as.character(result$ALIAS), split = "|", fixed = TRUE)
# table(sapply(l,length))
# 
# l1 = unique(do.call("c", l))
# 
# Rcpp:CharacterVector result = Rcpp:CharacterVector(l1.lenth()); 
# for(int i = 0; i < l1.length(); ++i)
# {
#     std::vector<std::string> res_i; 
#     for(int j = 0; j < l.length(); ++j)
#     {
#         for(int k = 0; k < l[j].length(); ++k)
#         {
#             if(l[j][k] == l1[i])
#             {
#                 res_i.push_back(j);
#             }
#         }
#     }
#     
#     if(res_i.length() == 0)
#     {
#         
#     } else {
#         
#     }
#     
# }
# length(l)
# length(1)
# 
# l2 = lapply(l1, function(x) lapply(l, function(y) x %in% y))
