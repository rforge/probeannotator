#' @keywords internal
#' @author Alexandre Thiery
#' @examples
#' ## using 'numeric' vector of length 3
#' vec = c(1,2,0.5)
#' ## this will work
#' .checkVector(vec, name = "vec", function.is="is.numeric", length.equal = 3, default = c(0,0,0))
#' f1 = function(vec) 
#'      .checkVector(vec, name = "vec", function.is="is.numeric", length.equal = 3, default = c(0,0,0))
#' f1(vec)
#' f1()
#' ## this will give and error
#' \dontrun{
#' ## * wrong length:
#' .checkVector(vec, name = "vec", function.is="is.numeric", length.equal = 4)
#' ## * wrong type:
#' .checkVector(vec, name = "vec", function.is="is.integer", length.equal = 3)
#' .checkVector(vec, name = "vec", function.is="is.character", length.equal = 3)
#' ## * 'vec' missing in environment:
#' .checkVector(   , name = "vec", function.is="is.numeric", length.equal = 3)
#' }
.checkVector <- function(vector, name, function.is, default, 
                         length.equal = 1, length.function = length,
                         name.is = sub("^is\\.","",function.is))
{
    if(missing(vector)) {
        if(missing(default)) 
        {
            stop(sprintf("'%s' is missing.",sQuote(name)))      
        }
        return(default)
    } else {
        
        ff = switch(function.is,
               is.integer=is.numeric,
               is.file=as.character,
               getFunction(function.is))
        
        if(ff(vector) == FALSE | length.function(vector)[1] != length.equal[1] )
        {
            stop( sprintf("%s must be an %s vector of length %d.",sQuote(name), name.is, length.equal) )  
        }         
        else 
        {
            if(function.is ==  "is.integer")
            {
                if(any(round(vector,0) != vector)) 
                    stop( sprintf("%s must be an %s vector of length %d.",sQuote(name), name.is, length.equal) )
            } else if(function.is ==  "is.file") 
            {
                if(!file.exists) 
                    stop( sprintf("%s doesn't exists.",sQuote(name)) )
            }
            return(vector)
        }
    }
}
