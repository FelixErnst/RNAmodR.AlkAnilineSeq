#' @include RNAmodR.AlkAnilineSeq.R
NULL

#' @name ModAlkAnilineSeq
#' @aliases AlkAnilineSeq ModAlkAnilineSeq ModSetAlkAnilineSeq
#' 
#' @author Felix G.M. Ernst [aut]
#' 
#' @title ModAlkAnilineSeq class to analyze AlkAnilineSeq data
#' 
#' @description 
#' 7-methyl guanosine (m7G), 3-methyl cytidine (m3C) and Dihydrouridine (D) 
#' are commonly found in rRNA and tRNA and can be detected classically by
#' primer extension analysis. However, since the modifications do not interfere
#' with Watson-Crick base pairing, a specific chemical treatment is employed
#' to cause strand breaks specifically at the modified positions.
#' 
#' This classical protocol was converted to a high throughput sequencing 
#' method call AlkAnilineSeq and allows modified position be detected by an
#' accumulation of 5'-ends at the N+1 position. Since the identify of the 
#' unmodified nucleotide is different for the three modified nucleotides, they
#' modification can be detected at the same time from the same samples.
#' 
#' The \code{ModAlkAnilineSeq} class uses the  
#' \code{\link[RNAmodR:NormEndSequenceData-class]{NormEnd5SequenceData}}
#' class to store and aggregate data along the transcripts. This includes 
#' normalized values against the whole transcript (normalzed cleavage) and 
#' normalized values against the overlapping reads (stop ratio), which are used
#' to score for modified positions.
#' 
#' In addition the 
#' \code{\link[RNAmodR:PileupSequenceData-class]{PileupSequenceData}} class is
#' used as well, to check, whether the base is can be called according to the
#' expected sequence identity.
#' 
#' Only samples named \code{treated} are used for this analysis. Normalization 
#' to untreated samples is currently not used.
#' 
#' @param x the input which can be of the different types depending on whether
#' a \code{ModRiboMethSeq} or a \code{ModSetRiboMethSeq} object is to be 
#' constructed. For more information have a look at the documentation of
#' the \code{\link[RNAmodR:Modifier-class]{Modifier}} and 
#' \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} classes.
#' @param annotation annotation data, which must match the information contained
#' in the BAM files. This is parameter is only required if \code{x} if not a 
#' \code{Modifier} object.
#' @param sequences sequences matching the target sequences the reads were 
#' mapped onto. This must match the information contained in the BAM files. This
#' is parameter is only required if \code{x} if not a \code{Modifier} object.
#' @param seqinfo An optional \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} 
#' argument or character vector, which can be coerced to one, to subset the 
#' sequences to be analyzed on a per chromosome basis.
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{minLength:} {The minimal read length to be used for the analysis 
#' (default: \code{minLength = 9L}).}
#' \item{minSignal:} {The minimal signal at the position as integer value 
#' (default: \code{minSignal = 10L}). If the reaction is very specific a lower
#' value may need to be used}
#' \item{minScoreNC:} {minimum for score (normalized cleavage) to identify m7G,
#' m3C and D positions de novo (default: \code{minScoreNC = 50L})}
#' \item{minScoreSR:} {minimum for score (stop ration) to identify m7G, m3C and D
#' positions de novo (default: \code{minScoreSR = 0.5})}
#' \item{minScoreBaseScore:} {minimum score for base calling (0.0-1.0)  
#' (default: \code{minScoreSR = 0.9})}
#' \item{scoreOperator:} {how the minimal score should be used as logical 
#' operator. "&" requires all minimal values to be exceeded, whereas "|" detects
#' positions, if at least one minimal values is exceeded (default: 
#' \code{scoreOperator = "&"}).}
#' \item{other arguments} {which are passed on to 
#' \code{\link[RNAmodR:EndSequenceData-class]{End5SequenceData}}}
#' }
#' 
#' @references 
#' - Marchand V, Ayadi L, __Ernst FGM__, Hertler J, Bourguignon-Igel V,
#' Galvanin A, Kotter A, Helm M, __Lafontaine DLJ__, Motorin Y (2018): 
#' "AlkAniline-Seq: Profiling of m7 G and m3 C RNA Modifications at Single 
#' Nucleotide Resolution." Angewandte Chemie (International ed. in English) 57 
#' (51), P. 16785–16790. DOI: 
#' \href{https://doi.org/10.1002/anie.201810946}{10.1002/anie.201810946}.
NULL

#' @rdname ModAlkAnilineSeq
#' @export
setClass("ModAlkAnilineSeq",
         contains = c("Modifier"),
         prototype = list(mod = c("m7G","m3C","D"),
                          score = "scoreNC",
                          dataType = c("NormEnd5SequenceData",
                                       "PileupSequenceData")))

# constructors -----------------------------------------------------------------

#' @rdname ModAlkAnilineSeq
#' @export
ModAlkAnilineSeq <- function(x, annotation = NA, sequences = NA, seqinfo = NA,
                             ...){
  RNAmodR::Modifier("ModAlkAnilineSeq", x = x, annotation = annotation, 
                    sequences = sequences, seqinfo = seqinfo, ...)
}

# settings ---------------------------------------------------------------------

#' @name ModAlkAnilineSeq-functions
#' @aliases aggregate modify settings visualizeData visualizeDataByCoord
#' 
#' @title Functions for ModAlkAnilineSeq
#' 
#' @description
#' All of the functions of \code{\link[RNAmodR:Modifier-class]{Modifier}} and
#' the \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} classes are
#' inherited by the \code{ModAlkAnilineSeq} and \code{ModSetAlkAnilineSeq}
#' classes.
#' 
#' @param x a \code{\link[RNAmodR:Modifier-class]{Modifier}} or a
#' \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} object. For more
#' details see also the man pages for the functions mentioned below.
#' @param value See \code{\link[RNAmodR:Modifier-class]{settings}}
#' @param force See \code{\link[RNAmodR:aggregate]{aggregate}}
#' @param coord,name,from,to,type,window.size,... See 
#' \code{\link[RNAmodR:visualizeData]{visualizeData}}.
#' 
#' @details 
#' \code{ModAlkAnilineSeq} specific arguments for \link{visualizeData}:
#' \itemize{
#' \item{\code{colour} - }{a named character vector of \code{length = 4} 
#' for the colours of the individual histograms. The names are expected to be 
#' \code{c("scoreNC","scoreSR")}}
#' }
#' 
#' @importMethodsFrom RNAmodR modify aggregate settings visualizeData 
#' visualizeDataByCoord
NULL

.norm_aas_args <- function(input){
  minLength <- 9L # for all scores
  minSignal <- 10L # for all scores
  minScoreNC <- 50L # for score normalized cleavage
  minScoreSR <- 0.50 # for score stop ratio
  minScoreBaseScore <- 0.9 # for base scoring
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
  if(!is.null(input[["minScoreSR"]])){
    minScoreSR <- input[["minScoreSR"]]
    if(!is.numeric(minScoreSR) | minScoreSR < 0 | minScoreSR > 1 ){
      stop("'minScoreSR' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreBaseScore"]])){
    minScoreBaseScore <- input[["minScoreBaseScore"]]
    if(!is.numeric(minScoreBaseScore) | minScoreBaseScore < 0 | 
       minScoreBaseScore > 1 ){
      stop("'minScoreBaseScore' must be numeric with a value between 0 and 1.",
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
            list(minLength = minLength,
                 minSignal = minSignal,
                 minScoreNC = minScoreNC,
                 minScoreSR = minScoreSR,
                 minScoreBaseScore = minScoreBaseScore,
                 scoreOperator = scoreOperator))
  args
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModAlkAnilineSeq"),
                 definition = function(x, value){
                   x <- callNextMethod()
                   value <- .norm_aas_args(value)
                   x@arguments[names(value)] <- unname(value)
                   x
                 })

# functions --------------------------------------------------------------------

.selected_base_score <- function(x,seq){
  seq <- IRanges::CharacterList(strsplit(as.character(seq),""))
  data <- x@unlistData
  # the order of baseSelect and baseColumns must be compatible
  # "conversion" T to U is performed with this
  baseSelect <- S4Vectors::DataFrame("G" = unlist(seq) == "G",
                                     "A" = unlist(seq) == "A",
                                     "U" = unlist(seq) == "U",
                                     "C" = unlist(seq) == "C")
  baseColumns <- which(colnames(data) %in% c("means.treated.G",
                                             "means.treated.A",
                                             "means.treated.T",
                                             "means.treated.C"))
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
  mod <- aggregate(sequenceData(x), condition = "Treated")
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
  score <- .selected_base_score(pileupData,sequences(x))
  # subset to columns needed
  data <- endData@unlistData[,c("means.treated.ends",
                                "means.treated.tx",
                                "means.treated.ol")]
  colnames(data) <- c("ends","scoreNC","scoreSR")
  # add base score
  data$baseScore <- score
  # split and return
  data <- IRanges::SplitDataFrameList(data)
  data@partitioning <- endData@partitioning
  data
}

#' @rdname ModAlkAnilineSeq-functions
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
      }
      x <- callNextMethod()
      x
    }
)

.get_aas_scores <- function(data){
  list(score = data$scoreNC,
       scoreSR = data$scoreSR)
}

.find_aas <- function(x){
  if(!hasAggregateData(x)){
    stop("Something went wrong.")
  }
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  grl <- ranges(x)
  # get the aggregate data
  mod <- aggregateData(x)
  # set up some arguments
  minSignal <- settings(x,"minSignal")
  minScoreNC <- settings(x,"minScoreNC")
  minScoreSR <- settings(x,"minScoreSR")
  minScoreBaseScore <- settings(x,"minScoreBaseScore")
  scoreOperator <- settings(x,"scoreOperator")
  # construct logical vector for passing the minSignal threshold
  signal <- IRanges::LogicalList(unname(mod@unlistData$ends >= minSignal))
  signal@partitioning <- mod@partitioning
  # find modifications
  modifications <- mapply(
    function(m,l,r,s){
      m <- m[!is.na(m$scoreNC) &
               !is.na(m$scoreSR) &
               !is.na(m$ends),,drop=FALSE]
      if(nrow(m) == 0L) return(NULL)
      m <- m[s,,drop=FALSE]
      if(nrow(m) == 0L) return(NULL)
      m <- m[mapply(Reduce,
                    rep(scoreOperator,nrow(m)),
                    m$scoreNC >= minScoreNC,
                    m$scoreSR >= minScoreSR,
                    m$baseScore >= minScoreBaseScore),,drop=FALSE]
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
    grl,
    signal)
  f <- !vapply(modifications,
               is.null,
               logical(1))
  modifications <- mapply(
    function(m,name){
      m$Parent <- rep(name,length(m))
      m
    },
    modifications[f],
    names(grl)[f],
    SIMPLIFY = FALSE)
  modifications <- GenomicRanges::GRangesList(modifications)
  message("done.")
  unname(unlist(modifications))
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod("modify",
          signature = c(x = "ModAlkAnilineSeq"),
          function(x, force = FALSE){
            if(!x@aggregateValidForCurrentArguments){
              x <- aggregate(x, force = TRUE)
            }
            x@modifications <- .find_aas(x)
            x <- callNextMethod()
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
                                seqinfo = NA, ...){
  RNAmodR::ModifierSet("ModAlkAnilineSeq", x, annotation = annotation,
                       sequences = sequences, seqinfo = seqinfo, ...)
}
