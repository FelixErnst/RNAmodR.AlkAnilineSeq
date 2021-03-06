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
#' \code{dataType} is \code{c("NormEnd3SequenceData","PileupSequenceData")}:
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
#' @return a \code{ModAlkAnilineSeq} or \code{ModSetAlkAnilineSeq} object
#' 
#' @references 
#' - Marchand V, Ayadi L, __Ernst FGM__, Hertler J, Bourguignon-Igel V,
#' Galvanin A, Kotter A, Helm M, __Lafontaine DLJ__, Motorin Y (2018): 
#' "AlkAniline-Seq: Profiling of m7 G and m3 C RNA Modifications at Single 
#' Nucleotide Resolution." Angewandte Chemie (International ed. in English) 57 
#' (51), P. 16785–16790. DOI: 
#' \href{https://doi.org/10.1002/anie.201810946}{10.1002/anie.201810946}.
#' 
#' @examples
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.AAS.gff3())
#' sequences <- RNAmodR.Data.example.AAS.fasta()
#' files <- list("wt" = c(treated = RNAmodR.Data.example.wt.1()),
#'               "Bud23del" = c(treated = RNAmodR.Data.example.bud23.1()),
#'               "Trm8del" = c(treated = RNAmodR.Data.example.trm8.1()))
#' # Creating a Modifier object of type ModRiboMethSeq
#' maas <- ModAlkAnilineSeq(files[[1]], annotation = annotation,
#'                          sequences = sequences)
#' # Creating a ModifierSet object of type ModSetRiboMethSeq
#' msaas <- ModSetAlkAnilineSeq(files, annotation = annotation,
#'                              sequences = sequences)
NULL

#' @rdname ModAlkAnilineSeq
#' @export
setClass("ModAlkAnilineSeq",
         contains = c("RNAModifier"),
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
#' @aliases aggregate modify settings plotData plotDataByCoord
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
#' @param coord,name,from,to,type,window.size,... See 
#' \code{\link[RNAmodR:plotData]{plotData}}.
#' 
#' @details 
#' \code{ModAlkAnilineSeq} specific arguments for \link{plotData}:
#' \itemize{
#' \item{\code{colour} - }{a named character vector of \code{length = 4} 
#' for the colours of the individual histograms. The names are expected to be 
#' \code{c("scoreNC","scoreSR")}}
#' }
#' 
#' @return 
#' \itemize{
#' \item{\code{settings}} {See 
#' \code{\link[RNAmodR:Modifier-functions]{settings}}.}
#' \item{\code{aggregate}} {See \code{\link[RNAmodR:aggregate]{aggregate}}.}
#' \item{\code{modify}} {See \code{\link[RNAmodR:modify]{modify}}.}
#' \item{\code{getDataTrack}} {a list of 
#' \code{\link[Gviz:DataTrack-class]{DataTrack}} object.}
#' \item{\code{plotData}} {See 
#' \code{\link[RNAmodR:plotData]{plotDataByCoord}}.}
#' \item{\code{plotDataByCoord}} {See 
#' \code{\link[RNAmodR:plotData]{plotDataByCoord}}.}
#' }
#' 
#' @importMethodsFrom RNAmodR modify aggregate settings plotData 
#' plotDataByCoord
#' 
#' @examples 
#' data(msaas,package="RNAmodR.AlkAnilineSeq")
#' maas <- msaas[[1]]
#' settings(maas)
#' aggregate(maas)
#' modify(maas)
#' getDataTrack(maas, "1", mainScore(maas))
NULL

.not_integer_bigger_than_zero <- RNAmodR:::.not_integer_bigger_than_zero
.not_numeric_bigger_zero <- RNAmodR:::.not_numeric_bigger_zero
.not_numeric_between_0_1 <- RNAmodR:::.not_numeric_between_0_1
.not_logical_operator <- RNAmodR:::.not_logical_operator

.ModAlkAnilineSeq_settings <- data.frame(
  variable = c("minLength",
               "minSignal",
               "minScoreNC",
               "minScoreSR",
               "minScoreBaseScore",
               "scoreOperator"),
  testFUN = c(".not_integer_bigger_than_zero",
              ".not_integer_bigger_than_zero",
              ".not_numeric_bigger_zero",
              ".not_numeric_between_0_1",
              ".not_numeric_between_0_1",
              ".not_logical_operator"),
  errorValue = c(TRUE,
                 TRUE,
                 TRUE,
                 TRUE,
                 TRUE,
                 TRUE),
  errorMessage = c("'minQuality' must be integer with a value higher than 0.",
                   "'minSignal' must be integer with a value higher than 0.",
                   "'minScoreNC' must be numeric with a value higher than 0.",
                   "'minScoreSR' must be numeric with a value between 0 and 1.",
                   "'minScoreBaseScore' must be numeric with a value between 0 and 1.",
                   "'scoreOperator' must be either '|' or '&'."),
  stringsAsFactors = FALSE)

.norm_aas_args <- function(input){
  minLength <- 9L # for all scores
  minSignal <- 10L # for all scores
  minScoreNC <- 50L # for score normalized cleavage
  minScoreSR <- 0.50 # for score stop ratio
  minScoreBaseScore <- 0.9 # for base scoring
  scoreOperator <- "&"
  args <- RNAmodR:::.norm_settings(input, .ModAlkAnilineSeq_settings, minLength,
                                   minSignal, minScoreNC, minScoreSR,
                                   minScoreBaseScore, scoreOperator)
  args <- c(RNAmodR:::.norm_Modifier_settings(input),
            args)
  args
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setReplaceMethod(f = "settings", 
    signature = signature(x = "ModAlkAnilineSeq"),
    definition = function(x, value){
        x <- callNextMethod()
        names <- names(value)
        value <- .norm_aas_args(value)
        x <- RNAmodR:::.add_settings_value(x, value, names)
        x
    }
)

# functions --------------------------------------------------------------------

.selected_base_score <- function(data,seq){
  seq <- IRanges::CharacterList(strsplit(as.character(seq),""))
  data <- data@unlistData
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

.apply_Np1_offset <- function(x){
  offsetAdd <- S4Vectors::DataFrame(as.list(rep(0,ncol(x@unlistData))))
  colnames(offsetAdd) <- colnames(unlist(x))
  x_rownames <- rownames(x)
  x <- x[IRanges::IntegerList(Map(seq.int,2L,lengths(x)))]
  x_offsetAdd <- IRanges::SplitDataFrameList(rep(list(offsetAdd),length(x)))
  x <- rbind(x,x_offsetAdd)
  rownames(x) <- x_rownames
  x
}

.aggregate_aas <- function(x){
  mod <- aggregate(sequenceData(x), condition = "Treated")
  if(is.null(names(mod))){
    names(mod) <- vapply(mod,class,character(1))
  }
  # get the NormEnd5SequenceData
  # shift each result one position upstream, since the signal has an N+1 offset
  endData <- .apply_Np1_offset(mod[["NormEnd5SequenceData"]])
  # get the PileupSequenceData
  # calculate base score
  pileupData <- mod[["PileupSequenceData"]]
  score <- .selected_base_score(pileupData,sequences(x))
  # subset to columns needed
  partitioning <- IRanges::PartitioningByEnd(endData)
  endData <- unlist(endData)[,c("means.treated.ends",
                             "means.treated.tx",
                             "means.treated.ol")]
  colnames(endData) <- c("ends","scoreNC","scoreSR")
  # add base score
  endData$baseScore <- score
  # split and return
  ans <- relist(endData, partitioning)
  rownames(ans) <- IRanges::CharacterList(RNAmodR:::.seqs_rl_strand(ranges(x)))
  ans
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "aggregateData", 
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = 
    function(x){
      .aggregate_aas(x)
    }
)

.get_aas_scores <- function(data){
  list(score = data$scoreNC,
       scoreSR = data$scoreSR)
}

.find_aas_modifications_in_data <- function(mod, letters, grl, signal,
                                            minScoreNC, minScoreSR,
                                            minScoreBaseScore,
                                            scoreOperator){
  modifications <- mapply(
    function(m,l,r,s){
      pos <- rownames(m)
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
      l <- l[which(pos %in% rownames(m))]
      if(!any(l %in% c("G","C","U"))){
        return(NULL)
      }
      ansm7G <- RNAmodR:::constructModRanges(
        r,
        m[l == "G",],
        modType = "m7G",
        scoreFun = RNAmodR.AlkAnilineSeq:::.get_aas_scores,
        source = "RNAmodR.AlkAnilineSeq",
        type = "RNAMOD")
      ansD <- RNAmodR:::constructModRanges(
        r,
        m[l == "U",],
        modType = "D",
        scoreFun = RNAmodR.AlkAnilineSeq:::.get_aas_scores,
        source = "RNAmodR.AlkAnilineSeq",
        type = "RNAMOD")
      ansm3C <- RNAmodR:::constructModRanges(
        r,
        m[l == "C",],
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
  f <- !vapply(modifications, is.null, logical(1))
  modifications <- mapply(
    function(m,name){
      m$Parent <- rep(name,length(m))
      m
    },
    modifications[f],
    names(grl)[f],
    SIMPLIFY = FALSE)
  modifications <- GenomicRanges::GRangesList(modifications)
  modifications
}

.find_aas_modifications <- function(x){
  if(!hasAggregateData(x)){
    stop("Something went wrong.")
  }
  grl <- ranges(x)
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  # get the aggregate data
  mod <- getAggregateData(x)
  # set up some arguments
  minSignal <- settings(x,"minSignal")
  minScoreNC <- settings(x,"minScoreNC")
  minScoreSR <- settings(x,"minScoreSR")
  minScoreBaseScore <- settings(x,"minScoreBaseScore")
  scoreOperator <- settings(x,"scoreOperator")
  # construct logical vector for passing the minSignal threshold
  signal <- mod[,"ends"] >= minSignal
  # find modifications
  modifications <- .find_aas_modifications_in_data(mod, letters, grl, signal,
                                                   minScoreNC, minScoreSR,
                                                   minScoreBaseScore,
                                                   scoreOperator)
  unname(unlist(modifications))
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod("findMod",
          signature = c(x = "ModAlkAnilineSeq"),
          function(x){
            .find_aas_modifications(x)
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
