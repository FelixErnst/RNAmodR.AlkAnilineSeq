#' @include RNAmodR.AlkAnilineSeq.R
NULL

#' @name ModAlkAnilineSeq
#' @aliases AlkAnilineSeq ModAlkAnilineSeq
#' 
#' @author Felix G.M. Ernst \email{felix.gm.ernst@@outlook.com}
#' 
#' @title ModAlkAnilineSeq
#' @description 
#' title
#' 
#' @param annotation
#' @param sequences
#' @param seqinfo
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
NULL

#' @rdname ModAlkAnilineSeq
#' @export
setClass("ModAlkAnilineSeq",
         contains = c("Modifier"),
         prototype = list(mod = c("m7G","m3C","D"),
                          score = "scoreNC",
                          dataType = c("NormEnd5SequenceData",
                                       "PileupSequenceData")))

# settings ---------------------------------------------------------------------
.norm_aas_args <- function(input){
  maxLength <- NA # for all scores
  minLength <- 10L # for all scores
  minSignal <- 10L # for all scores
  minScoreNC <- 50L # for score normalized cleavage
  minScoreSR <- 0.50 # for score stop ratio
  scoreOperator <- "&"
  if(!is.null(input[["minLength"]])){
    minLength <- input[["minLength"]]
    if(!is.integer(minLength) | minLength < 1L){
      if(!is.na(minLength)){
        stop("'minQuality' must be integer with a value higher than 0.",
             call. = FALSE)
      }
    }
  }
  if(!is.null(input[["minSignal"]])){
    minSignal <- input[["minSignal"]]
    if(!is.integer(minSignal) | minSignal < 1L){
      stop("'minSignal' must be integer with a value higher than 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreNC"]])){
    minScoreNC <- input[["minScoreNC"]]
    if(!is.numeric(minScoreNC) | minScoreNC < 0 ){
      stop("'minScoreNC' must be numeric with a value higher than 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreSR"]])){
    minScoreSR <- input[["minScoreSR"]]
    if(!is.numeric(minScoreSR) | minScoreSR < 0 | minScoreSR > 1 ){
      stop("'minScoreSR' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["scoreOperator"]])){
    scoreOperator <- input[["scoreOperator"]]
    if(!(scoreOperator %in% c("|","&"))){
      stop("'scoreOperator' must be either '|' or '&'.",
           call. = FALSE)
    }
  }
  args <- RNAmodR:::.norm_args(input)
  args <- c(args,
            list(maxLength = maxLength,
                 minLength = minLength,
                 minSignal = minSignal,
                 minScoreNC = minScoreNC,
                 minScoreSR = minScoreSR,
                 scoreOperator = scoreOperator))
  args
}

#' @name ModAlkAnilineSeq
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModAlkAnilineSeq"),
                 definition = function(x, value){
                   x <- callNextMethod()
                   value <- .norm_aas_args(value)
                   x@arguments[names(value)] <- unname(value)
                   x
                 })

# constructors -----------------------------------------------------------------

#' @rdname ModAlkAnilineSeq
#' @export
ModInosine <- function(x, annotation = NA, sequences = NA, seqinfo = NA, ...){
  Modifier("ModAlkAnilineSeq", x = x, annotation = annotation, 
           sequences = sequences, seqinfo = seqinfo, ...)
}

# functions --------------------------------------------------------------------

.calculate_base_score <- function(x,seq){
  seq <- IRanges::CharacterList(strsplit(as.character(seq),""))
  data <- x@unlistData
  baseSelect <- S4Vectors::DataFrame("G" = unlist(seq) == "G",
                          "A" = unlist(seq) == "A",
                          "U" = unlist(seq) == "U",
                          "C" = unlist(seq) == "C")
  baseColumns <- which(colnames(data) %in% c("means.G",
                                             "means.A",
                                             "means.T",
                                             "means.C"))
  baseData <- as(data[,baseColumns],"NumericList")
  baseSelect <- as(baseSelect,"LogicalList")
  tmpBaseValues <- unname(baseData)[unname(baseSelect)]
  baseValues <- rep(0,nrow(data))
  baseValues[which(baseSelect[[1]])] <- tmpBaseValues[[1]]
  baseValues[which(baseSelect[[2]])] <- tmpBaseValues[[2]]
  baseValues[which(baseSelect[[3]])] <- tmpBaseValues[[3]]
  baseValues[which(baseSelect[[4]])] <- tmpBaseValues[[4]]
  baseValues[is.na(baseValues)] <- 0
  baseValues[is.infinite(baseValues)] <- 0
  baseValues
}

.aggregate_aas <- function(x){
  message("Aggregating data and calculating scores...")
  mod <- aggregate(seqData(x),
                   condition = "Treated")
  if(is.null(names(mod))){
    names(mod) <- vapply(mod,class,character(1))
  }
  # get the NormEnd5SequenceData
  # shift each result one position upstream, since the signal has an N+1 offset
  endData <- mod[["NormEnd5SequenceData"]]
  offsetAdd <- S4Vectors::DataFrame(as.list(rep(0,ncol(endData@unlistData))))
  colnames(offsetAdd) <- colnames(endData@unlistData)
  endData <- IRanges::SplitDataFrameList(lapply(endData,
                                                function(m){
                                                  m <- m[seq.int(2L,nrow(m)),]
                                                  rbind(m,offsetAdd)
                                                }))
  # get the PileupSequenceData
  # calculate base score
  pileupData <- mod[["PileupSequenceData"]]
  score <- .calculate_base_score(pileupData,sequences(x))
  # subset to columns needed
  data <- endData@unlistData[,c("means.ends",
                                "means.tx",
                                "means.ol")]
  colnames(data) <- c("ends","scoreNC","scoreSR")
  # add base score
  data$baseScore <- score
  # split and return
  data <- IRanges::SplitDataFrameList(data)
  data@partitioning <- endData@partitioning
  data
}

#' @name ModAlkAnilineSeq
#' @export
setMethod(
  f = "aggregate", 
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = 
    function(x, force = FALSE){
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

.get_aas_scores <- function(data){
  list(score = data$scoreNC,
       scoreSR = data$scoreSR)
}

.find_aas <- function(x){
  message("Searching for m7G/m3C/D ...")
  #
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  ranges <- ranges(x)
  ranges <- split(RNAmodR:::.get_parent_annotations(ranges),
                  seq_along(ranges))
  # get the aggregate data
  mod <- aggregateData(x)
  # set up some arguments
  minSignal <- settings(x,"minSignal")
  minScoreNC <- settings(x,"minScoreNC")
  minScoreSR <- settings(x,"minScoreSR")
  scoreOperator <- settings(x,"scoreOperator")
  # construct logical vector for passing the minSignal threshold
  signal <- IRanges::LogicalList(unname(mod@unlistData$ends >= minSignal))
  signal@partitioning <- mod@partitioning
  # find modifications
  modifications <- mapply(
    function(m,l,r,s){
      rownames(m) <- seq_len(BiocGenerics::width(r))
      m <- m[!is.na(m$scoreNC) &
               !is.na(m$scoreSR) &
               !is.na(m$ends),]
      if(nrow(m) == 0L) return(NULL)
      m <- m[s,]
      if(nrow(m) == 0L) return(NULL)
      m <- m[mapply(Reduce,
                    rep(scoreOperator,nrow(m)),
                    m$scoreNC >= minScoreNC,
                    m$scoreSR >= minScoreSR,
                    m$baseScore >= 0.9),]
      if(nrow(m) == 0L) return(NULL)
      ansm7G <- RNAmodR:::.constructModRanges(
        r,
        m[l[as.integer(rownames(m))] == "G",],
        modType = "m7G",
        scoreFun = RNAmodR.AlkAnilineSeq:::.get_aas_scores,
        source = "RNAmodR.AlkAnilineSeq",
        type = "RNAMOD")
      ansD <- RNAmodR:::.constructModRanges(
        r,
        m[l[as.integer(rownames(m))] == "U",],
        modType = "D",
        scoreFun = RNAmodR.AlkAnilineSeq:::.get_aas_scores,
        source = "RNAmodR.AlkAnilineSeq",
        type = "RNAMOD")
      ansm3C <- RNAmodR:::.constructModRanges(
        r,
        m[l[as.integer(rownames(m))] == "C",],
        modType = "m3C",
        scoreFun = RNAmodR.AlkAnilineSeq:::.get_aas_scores,
        source = "RNAmodR.AlkAnilineSeq",
        type = "RNAMOD")
      ans <-  unlist(GenomicRanges::GRangesList(list(ansm7G,ansD,ansm3C)))
      ans <- ans[order(BiocGenerics::start(ans))]
      ans
    },
    mod,
    letters,
    ranges,
    signal)
  modifications <- GenomicRanges::GRangesList(
    modifications[!vapply(modifications,
                          is.null,
                          logical(1))])
  unname(unlist(modifications))
}

#' @rdname ModAlkAnilineSeq
#' @export
setMethod("modify",
          signature = c(x = "ModAlkAnilineSeq"),
          function(x, force){
            # get the aggregate data
            x <- aggregate(x, force)
            x@modifications <- .find_aas(x)
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
ModSetAlkAnilineSeq <- function(x, annotation = NA, sequences = NA,
                                seqinfo = NA){
  ModifierSet("ModAlkAnilineSeq", x, annotation = annotation,
              sequences = sequences, seqinfo = seqinfo)
}
