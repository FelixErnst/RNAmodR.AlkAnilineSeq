#' @include RNAmodR.AlkAnilineSeq.R
NULL

#' @name ModAlkAnilineSeq
#' @aliases AlkAnilineSeq ModAlkAnilineSeq
#' 
#' @title ModAlkAnilineSeq
#' @description 
#' title
#' 
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{minSignal:}{The minimal singal at the position as integer value 
#' (default: \code{minSignal = 10L}). If the reaction is very specific a lower
#' value may need to be used}
#' \item{minScore:}{minimum for score to identify m7G, m3C and D positions 
#' de novo (default: \code{minScore = 0.5})}
#' \item{maxLength:}{The default read length. Reads with this length or longer
#' are discarded, since they represent non-fragemented reads. This might need to
#' be adjusted for individual samples dending on the experimental conditions.
#' This is argument is passed on to \code{\link{ProtectedEndSequenceData}} 
#' (default: \code{maxLength = 50L})}
#' \item{other arguments}{which are passed on to 
#' \code{\link{End5SequenceData}}}
#' }
#' 
#' 
NULL

#' @rdname ModAlkAnilineSeq
#' @export
setClass("ModAlkAnilineSeq",
         contains = c("Modifier"),
         prototype = list(mod = c("m7G","m3C","D"),
                          score = "score",
                          dataClass = "End5SequenceData"))


setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModAlkAnilineSeq"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              fasta,
                              gff)
    return(.Object)
  }
)

# constructors -----------------------------------------------------------------

.norm_aas_args <- function(input){
  maxLength <- 50L # for all scores
  minSignal <- 10L # for all scores
  minScore <- 0.75 # for score C/RMS
  if(!is.null(input[["minSignal"]])){
    minSignal <- input[["minSignal"]]
    if(!is.integer(minSignal) | minSignal < 1L){
      stop("'minSignal' must be integer with a value higher than 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScore"]])){
    minScore <- input[["minScore"]]
    if(!is.numeric(minScore) | minScore < 0 | minScore > 1){
      stop("'minScore' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  args <- .norm_args(input)
  args <- c(args,
            list(maxLength = maxLength,
                 minSignal = minSignal,
                 minScore = minScore))
  args
}

#' @name ModAlkAnilineSeq
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModAlkAnilineSeq"),
                 definition = function(x,value){
                   x <- callNextMethod()
                   value <- .norm_aas_args(value)
                   x@arguments[names(value)] <- unname(value)
                   x
                 })


setGeneric( 
  name = "ModAlkAnilineSeq",
  def = function(x,
                 ...) standardGeneric("ModAlkAnilineSeq")
)
# Create Modifier class from file character, fasta and gff file
#' @rdname ModAlkAnilineSeq
#' @export
setMethod("ModAlkAnilineSeq",
          signature = c(x = "character"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            ans <- RNAmodR:::.ModFromCharacter("ModAlkAnilineSeq",
                                               x,
                                               fasta,
                                               gff,
                                               list(...))
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 list(...))
            ans
          }
)

# Create Modifier class from bamfiles, fasta and gff file
#' @rdname ModAlkAnilineSeq
#' @export
setMethod("ModAlkAnilineSeq",
          signature = c(x = "BamFileList"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            ans <- RNAmodR:::.ModFromCharacter("ModAlkAnilineSeq",
                                               x,
                                               fasta,
                                               gff,
                                               list(...))
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 list(...))
            ans
          }
)

# Create Modifier class from existing SequenceData
#' @rdname ModAlkAnilineSeq
#' @export
setMethod("ModAlkAnilineSeq",
          signature = c(x = "SequenceData"),
          function(x,
                   modifications = NULL,
                   ...){
            ans <- RNAmodR:::.ModFromSequenceData("ModAlkAnilineSeq",
                                                  x,
                                                  list(...))
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 list(...))
            ans
          }
)

# functions --------------------------------------------------------------------

# calculates score C according to Birkedal et al. 2015
# it is simplified to use the weighted neighboring positions
# and not the left and right area seperatly
.calculate_alkaniline_score_c <- function(n,
                                             areaL,
                                             weightsL,
                                             areaR,
                                             weightsR){
  waL <- sum(weightsL * areaL) / sum(weightsL)
  waR <- sum(weightsR * areaR) / sum(weightsR)
  dividend <- n
  divisor <- 0.5 * ( waL + waR )
  ans <- 1 - (dividend / divisor)
  ans <- vapply(ans,max,numeric(1),0)
  ans[is.na(ans)] <- 0
  return(ans)
}
.calculate_alkaniline_score <- compiler::cmpfun(.calculate_alkaniline_score_c)

# AlkAnilineSeq scores ---------------------------------------------------------

.get_score_alkanilineseq <- function(data,
                            countsL,
                            weightsL,
                            countsR,
                            weightsR){
  scoreMeth <- NumericList(mapply(FUN = .calculate_alkaniline_score,
                                  data,
                                  countsL,
                                  weightsL,
                                  countsR,
                                  weightsR))
  return(scoreMeth)
}

.aggregate_aas <- function(x){
  browser()
  message("Aggregating data and calculating scores...")
  # get the means. the sds arecurrently disregarded for this analysis
  mod <- aggregate(seqData(x),
                   condition = "Treated")
  means <- IntegerList(mod@unlistData[,which(grepl("mean",
                                                   colnames(mod@unlistData)))])
  means@partitioning <- mod@partitioning
  # set up variables
  n <- length(mod)
  nV <- seq_len(n)
  lengths <- lengths(mod)
  pos <- lapply(lengths,seq_len)
  
  
  
  # calculate the actual scores
  score <- .get_score_alkanilineseq(means,
                                    neighborCountsL,
                                    weightsListL,
                                    neighborCountsR,
                                    weightsListR)
  # scoreMAX <- .get_score_max()
  ans <- DataFrame(ends = unlist(means),
                   score = unlist(score),
                   row.names = NULL)
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- mod@partitioning
  ans
}

#' @name ModAlkAnilineSeq
#' @export
setMethod(
  f = "aggregate", 
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = 
    function(x,
             force = FALSE){
      if(missing(force)){
        force <- FALSE
      }
      if(!hasAggregateData(x) || force){
        x@aggregate <- .aggregate_aas(x)
        x@aggregateValidForCurrentArguments <- TRUE
      }
      x
    }
)

.get_rms_scores <- function(data){
  list(score = data$scoreRMS,
       scoreA = data$scoreA,
       scoreB = data$scoreB)
}

.find_rms <- function(x){
  message("Searching for m7G/m3C/D ...")
  #
  data <- seqData(x)
  letters <- CharacterList(strsplit(as.character(sequences(data)),""))
  ranges <- split(.get_parent_annotations(ranges(data)),
                  seq_along(ranges(data)))
  # get the aggregate data
  mod <- aggregateData(x)
  # setup args
  minSignal <- settings(x,"minSignal")
  minScoreRMS <- settings(x,"minScore")
  # find modifications
  modifications <- mapply(
    function(m,l,r){
      rownames(m) <- seq_len(width(r)) + 1
      m <- m[!is.na(m$scoreA) &
               !is.na(m$scoreB) &
               !is.na(m$scoreRMS),]
      if(nrow(m) == 0L) return(NULL)
      m <- m[mapply(Reduce,
                    rep(scoreOperator,nrow(m)),
                    m$scoreA >= minScoreA,
                    m$scoreB >= minScoreB,
                    m$scoreRMS >= minScoreRMS),]
      m <- m[m$ends >= minSignal,]
      if(nrow(m) == 0L) return(NULL)
      ans <- .construct_mod_ranges(r,
                                   m,
                                   modType = "Am",
                                   scoreFun = .get_rms_scores,
                                   source = "RNAmodR",
                                   type = "RNAMOD")
      ans$mod <- paste0(l[start(ans)],"m")
      ans
    },
    mod,
    letters,
    ranges)
  modifications <- GRangesList(modifications[!vapply(modifications,
                                                     is.null,
                                                     logical(1))])
  unname(unlist(modifications))
}

#' @rdname ModAlkAnilineSeq
#' @export
setMethod("modify",
          signature = c(x = "ModAlkAnilineSeq"),
          function(x,
                   force){
            # get the aggregate data
            x <- aggregate(x, force)
            x@modifications <- .find_rms(x)
            x@modificationsValidForCurrentArguments <- TRUE
            message("done.")
            x
          }
)

# ModSetAlkAnilineSeq ----------------------------------------------------------

#' @rdname ModAlkAnilineSeq
#' @export
setClass("ModSetAlkAnilineSeq",
         contains = "ModifierSet",
         prototype = list(elementType = "ModAlkAnilineSeq"))

#' @rdname ModAlkAnilineSeq
#' @export
ModSetAlkAnilineSeq <- function(x, fasta = NA, gff = NA){
  ModifierSet("ModAlkAnilineSeq", x, fasta = fasta, gff = gff)
}

