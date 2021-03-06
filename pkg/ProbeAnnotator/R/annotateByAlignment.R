#' Fast Gene Annotation by Alignment or Chromosomal Location
#' @rdname annotateBy
#' @description These method allows to perform fast probe-to-gene annotation using \emph{1)} chromosomal location, or \emph{2)} alignment files from softwares such \code{"bowtie"}, \code{"bowtie2"} or \code{"gmap"}.
#' 
#' @param x A \code{GRanges} object or \code{data.frame}, giving the coordinates to annotate.
#' @param file1 A character vector, the name of the probes's fasta file. See \emph{Input Files}.
#' @param probesetSep A character vector, indicating the probeset seperator if probes are organised in sets.
#' @param file2 A character vector, the name of the probes-to-reference SAM file. See \emph{Input Files}.
#' @param sepFile2 A character vector, the string column separator un \code{file2}.See \emph{Input Files}.
#' @param alignment.columnsIndex A numeric vector, containing the index of the \emph{score}, \emph{probe's name}, \emph{reference's name} and \emph{alignment offset} in the alignment file \code{file2}. This argument is not used if \code{alignment.method} is declared.
#' @param minScore A numeric value, giving the minimum allowed alignment score.
#' @param refUpStream A numeric value, giving the number of upstream bp in the reference, see \emph{Reference Sequence Region}.
#' @param refDownStream A numeric value, giving the number of downstream bp in the reference, see \emph{Reference Sequence Region}. 
#' @param mapType A character vector representing the probe-to-gene mapping type. This must be one of \code{"EXONINTRON"}, \code{"NO_EXONINTRON"} or \code{"EXON"}. Any unambiguous substring can be given. See \emph{Reference Sequence Region}.
#' @param txDb A \code{TranscriptDb} object, giving the genomic references for the alignement. If \code{txDb} is missing, the default \code{TxDb.Hsapiens.UCSC.hg19.knownGene} will be used.
#' @param orgDb A \code{OrganismDb} object, giving the details of the genomic references. If \code{orgDb} is missing, the default \code{org.Hs.eg.db} will be used.
#' @param orgDb_Columns A character vector (optional), giving which columns to extract from the \code{orgDb}. Note that if \code{orgDb_Columns} is used, then the user selected columns will be selected instead of the default request on \code{org.Hs.eg.db}, see details section.
#' @param promotorRange A integer vector, giving the window size for the genes' promotor site in bp. Default is \code{1500}, see \emph{Location Types}.
#' @param extendedRange A integer vector, giving the window size for the genes' extended site in bp. Default is \code{2000}, see \emph{Location Types}.
#' @param sep_intra A character vector, giving the separator character for gene information, see \emph{Ouput Separators}.
#' @param sep_inter A character vector, giving the separator character between genes, see \emph{Ouput Separators}.
#' @param verbose A logical value, indicating if messages should be printed. Default is \code{FALSE}.
#' 
#' @section Reference Sequence Region:
#' 
#' \describe{
#' \item{Genomic level:}{it corresponds to the gene reference format used for alignment and is controled with the \code{mapType} argument. It determines the \code{loctype} in the annotation. There are three types of allowed \code{mapType}:
#' \enumerate{
#' \item \code{"NO_EXONINTRON"}: if the reference sequence contains one sequence per gene. 
#' 
#' The available \code{loctype} are: \code{"gene"},\code{"promotor"},\code{"extended"},\code{"intragenic"}.
#' \item \code{"EXONINTRON"}: if the reference contains one sequence per transcript, including introns.
#'
#' The available \code{loctype} are: \code{"gene"},\code{"intron"},\code{"exon"},\code{"promotor"},\code{"extended"},\code{"intragenic"}.
#' \item \code{"EXON"}: if the reference contains one sequence per transcript, without introns.
#'  
#' The available \code{"gene"},\code{"exon"},\code{"promotor"},\code{"extended"},\code{"intragenic"}.}}
#' 
#' \item{Upstream and downstream.}{The upstream and downstream values that are retrieved in the reference are controled with the \code{refUpStream} and \code{refDownStream} parameters (in bp). }
#' }
#' 
#' @section Input Files:
#' 
#' 
#' \describe{
#' \item{Probes' FASTA file:}{this file's name is given by the \code{file1} argument. It is used to retrieve all of the platform's probe names.
#' }
#' #' \item{Probeset:}{when the probes in the platform are arranged in probesets, one can use the \code{probesetSep} to define the probesets seperator string. 
#' 
#' For example, using Affymetrix's XXX platform, set \code{probesetSep="at."}.
#' }
#' }
#' 
#' \describe{
#' \item{Alignement output:}{this file's name is given by the \code{file2} argument. Those outputs must have columns the alignment \emph{score}, \emph{probe's name}, \emph{reference's name} and \emph{alignment offset} (see \emph{Alignment format} below).
#'
#' }
#' \item{Alignment format:}{the alignment format must be known to this function to get the alignment infomration (score, probe, ref ,offset). The default input is the SAM format (see specifications at \url{https://samtools.github.io/hts-specs/SAMv1.pdf}, however it can be achieved manualy using the \code{alignment.columnsIndex} and the \code{sepFile2} arguments.
#' 
#' \code{alignment.columnsIndex} and the \code{sepFile2} allows user to enter specific alignment ouput format. The \code{alignment.columnsIndex} must indicate the columns of \emph{score}, \emph{probe's name}, \emph{reference's name} and \emph{alignment offset}. The column separator is given with argument \code{sepFile2}.}
#' }
#' @section Location Types:
#' 
#' The location types (column \code{loctype}) are pre-defined regions that describe gene's region to which the probe match to. There are six types or \code{loctype} (shown in table below).
#' 
#' \tabular{rc}{
#' \code{loctype} \tab illustration \cr
#' \code{"gene"} \tab     \code{.............########################.............} \cr
#' \code{"intron"} \tab   \code{................***.....******..**................} \cr
#' \code{"exon"} \tab     \code{.............===...=====......==..===.............} \cr
#' \code{"promotor"} \tab \code{..........<+++++>.................................}\cr
#' \code{"extended"} \tab \code{.......<~~~~>........................<~~~~>.......}\cr
#' \code{"intragenic"} \tab \code{------>....................................<------}
#' }
#' 
#' 
#' \describe{
#' \item{Promotor and extended regions:}{they can be adjusted using the \code{promotorRange} and \code{extendedRange} parameters. The \emph{promotor}'s range is set at \eqn{\pm}{+/-} \code{promotorRange} bp from the gene's start location. The \emph{extended}'s ranges are located at both ends of the gene, extending the gene region by \code{extendedRange} bp.
#' 
#' To exclude the \code{"promotor"} and/or \code{"extended"} regions for the annotation, set \code{promotorRange=0} and/or \code{extendedRange=0}.}
#' }
#' 
#' @section Ouput Separators:
#' \describe{
#' \item{\code{sep_intra}}{controls the columns' elements concatenation in a unique reference (\emph{i.e.} genes). If \code{sep_intra=";"} then gene items that have multiple entries are concatenated with \code{";"}. For example:
#' \itemize{
#' \item EGFR gene has six other symbols (ERBB, HER1, mENA, ERBB1, PIG61 and NISBD2), the \code{"alias"} column will be:
#' 
#' \code{EGFR;ERBB;HER1;mENA;ERBB1;PIG61;NISBD2}
#' \item HOXA10 gene has four other symbols (HOX1, HOX1.8, HOX1H and PL), the \code{"alias"} column will be:
#' 
#' \code{HOXA10;HOX1;HOX1.8;HOX1H;PL}
#' }} 
#' \item{\code{sep_inter}}{controls the columns' elements concatenation when a probe is mapped to multiple references (\emph{i.e.} genes). 
#' 
#' For example assume that a probe is mapped both to HOXA10 and EGFR, then all columns containing gene information are concatenated with \code{sep_inter="//"}. Here:
#' 
#' \code{symbol   alias}
#' 
#' \code{HOXA10//EGFR   HOXA10;HOX1;HOX1.8;HOX1H;PL//EGFR;ERBB;HER1;mENA;ERBB1;PIG61;NISBD2}
#' }}
#' @return
#' The method \code{annotateByAlignment} returns a \code{data.frame} object, with one record (or row) for each probes given in \code{x}. 
#' 
#' With the default organism database, this \code{data.frame} contains the following information:
#' \tabular{rl}{
#' \strong{Column} \tab \strong{Comment} \cr
#' \code{probe_name} \tab xxxx \cr
#' \code{entrezid} \tab xxxx \cr
#' \code{chr} \tab xxxx \cr
#' \code{strand} \tab xxxx \cr
#' \code{loctype} \tab xxxx \cr
#' \code{gene_end} \tab xxxx \cr
#' \code{gene_start} \tab xxxx \cr
#' \code{gene_symbol} \tab xxxx \cr
#' \code{gene_alias} \tab xxxx \cr 
#' \code{gene_name} \tab xxxx
#' }
#' 
#' If the user supplies its own organism database (\code{orgDb} and \code{orgDb_Columns}), the function will return a equivalent \code{data.frame} as above, with the columns \code{gene_symbol}, \code{gene_symbol} and \code{gene_symbol} replaced by the ones provided in \code{orgDb_Columns}.
#' @examples
#' ## .. todo
#' @export
#' @import GenomicFeatures
#' @importFrom GenomicRanges as.data.frame
#' @importFrom AnnotationDbi select
#' @importFrom RSQLite dbGetQuery
#' @import Rcpp
#' @useDynLib ProbeAnnotator
annotateByAlignment = function(file1, file2,
                              alignment.columnsIndex, sepFile2, 
                              minScore,         
                              refDownStream, refUpStream,
                              probesetSep,
                              txDb, 
                              mapType = "EXONINTRON", 
                              promotorRange = 1500, extendedRange = 2000, 
                              orgDb, orgDb_Columns, 
                              sep_intra = ";", sep_inter = '//', verbose = FALSE)
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
    
    ## check refDownStream & refUpStream
    refDownStream = .checkVector(refDownStream, "is.integer", name = quote(refDownStream))
    refUpStream = .checkVector(refUpStream, "is.integer", name = quote(refUpStream))
    ## check file1, file2, sepFile2, 
    file1 = .checkVector(file1, "is.file", name = quote(file1))
    file2 = .checkVector(file2, "is.file", name = quote(file2))
    sepFile2 = .checkVector(sepFile2, "is.character", name = quote(sepFile2))
    ## check minScore
    minScore = .checkVector(minScore, "is.integer", name = quote(minScore), default = -1)
    
    ## alignment.method & alignment.columnsIndex, 
    if(!missing(alignment.method))
    {
        alignment.method = .checkVector(alignment.method, "is.character", name = quote(alignment.method))
        
        match1 = match.arg(alignment.method, choices = c("bowtie", "bowtie2","gmap"), several.ok = FALSE)
        alignment.columnsIndex = switch(match1, 
               bowtie=c(5,1,3,4)-1,
               bowtie2=c(5,1,3,4)-1,
               gmap=c(0,-1,-1,0),
               stop("Unknown alignment method."))
        
        sepFile2 = switch(match1, 
                          bowtie="\t",
                          bowtie2="\t",
                          gmap="\t",
                          stop("Unknown alignment method."))
    } else if(!missing(alignment.columnsIndex) & !missing(sepFile2))
    {
        ## check probesetSep
        probesetSep = .checkVector(probesetSep, "is.character", name = quote(probesetSep), default = "")
        ## check alignment.columnsIndex
        alignment.columnsIndex = .checkVector(alignment.columnsIndex, "is.integer", name = quote(alignment.columnsIndex), 
                                              length.equal = 4)
        w = alignment.columnsIndex[1:3] < 0
        if(any(w))
            sprintf("Column(s) index must be positive (%s < 0).", paste(c("score","probe","reference")[w], collapse = ","))
    } else {
        stop("'alignment.columnsIndex' and 'sepFile2' or 'alignment.method' must be specified.")
    }
   
    
    ## set x
    x = list(
        c(file1, file2, sepFile2,probesetSep),
        c(unlist(alignment.columnsIndex), minScore, refDownStream, refUpStream)
        )
    
    refList= annotateGetRefFromTranscriptDb(txDb, exons = exons, verbose = verbose)
    infList= annotateGetInfoFromOrganismtDb(orgDb, orgDb_Columns, sep_intra, verbose = verbose)
    
    minRange = max(c(promotorRange, extendedRange, 0))
    maxRange = max(c(extendedRange, 0))
    txAnnotRange = c(promotorRange, extendedRange, minRange, maxRange)
    txAnnotRange = as.integer(txAnnotRange)
    
    if(verbose) message("Perform annotation...")
    result <- .Call('ProbeAnnotator_annotateByAlignment', PACKAGE = 'ProbeAnnotator',
                    x, refList$df_GENETR, refList$txGroup, refList$df_EXON, refList$stack_EXON, 
                    infList$txInfo, infList$orgDb_Columns,
                    c(sep_intra,sep_inter),                
                    txAnnotRange, 
                    list( c(mapTypeIndex, as.integer(verbose)) )
                    )
    
    
    if(verbose) message("Done")
    if(!is.null(result))
        return(as.data.frame(do.call("cbind", result)))    
}