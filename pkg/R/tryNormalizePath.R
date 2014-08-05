#' @title tryNormalizePath
#'
#' A function check that files or directories exists, and normalize their path.
#'
#' @name tryNormalizePath
#' @param path a character vector giving the name directory or file to check.
#' @param is.file a logical value indicating wether the \code{path} is a directory or a file.
#' @param check.file a logical value indicating (in the case of a file, when \code{is.file} is \code{TRUE}), if the file should be checked (\code{check.file=TRUE}), or if the directory should be checked (\code{check.file=FALSE})
tryNormalizePath <- function(path, is.file = TRUE, check.file=TRUE) {
  
  is.file.str <- switch(as.character(is.file),
         "TRUE" = "file",
         "FALSE" = "directory")
  
  if(is.null(path)) {
    message(sprintf("Error: input '%s' is NULL", is.file.str))
    return(NULL)
  }
  
  result <- tryCatch({
    base <- NULL
    if(is.file & !check.file)
    {
      base <- basename(path)
      path <- dirname(path)
    }
    new.path <- normalizePath(path,mustWork=TRUE)

    if(is.file & !check.file)
    {
      new.path <- file.path(path, base)
    }
    return(new.path)
  }, error = function(error) {
    message(sprintf("Error: the %s '%s' does not exists.", is.file.str, path))
  })
  return(result)
}

