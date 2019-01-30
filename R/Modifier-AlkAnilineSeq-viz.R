#' @include RNAmodR.AlkAnilineSeq.R
#' @include Modifier-AlkAnilineSeq-class.R
NULL

RNAMODR_AAS_PLOT_DATA <- c("ends",
                           "scoreNC",
                           "scoreSR")
RNAMODR_AAS_PLOT_DATA_NAMES <- c(ends = "5'-ends",
                                 scoreNC = "Normalized Cleavage",
                                 scoreSR = "Stop ratio")
RNAMODR_AAS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   scoreNC = "#DECBE4",
                                   scoreSR = "#CCEBC5")
#' @rdname ModAlkAnilineSeq
#' 
#' @name visualizeData
#' @details 
#' \code{ModAlkAnilineSeq} specific arguments for \link{visualizeData}:
#' \itemize{
#' 
#' }
NULL

#' @rdname ModAlkAnilineSeq
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModAlkAnilineSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("ends","scoreNC","scoreSR"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA
    }
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModAlkAnilineSeq
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = function(x, name, from, to,
                        type = c("ends","scoreNC","scoreSR"), ...) {
    browser()
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA
    }
    callNextMethod(x = x, name, from,  to, type = type, ...)
  }
)

setMethod(
  f = ".dataTracks",
  signature = signature(x = "ModAlkAnilineSeq",
                        data = "GRanges",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, data, seqdata, sequence, args) {
    n <- ncol(mcols(data))
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_AAS_PLOT_DATA_COLOURS[seq.int(1,n)]
    }
    if(is.null(names(colour))){
      names(colour) <- colnames(mcols(data))
    }
    dts <- lapply(seq_len(n),
                  function(i){
                    column <- colnames(mcols(data)[i])
                    colour <- colour[column]
                    name <- RNAMODR_AAS_PLOT_DATA_NAMES[column]
                    dt <- Gviz::DataTrack(data,
                                          data = column,
                                          name = name,
                                          fill = colour,
                                          col = colour,
                                          type = "histogram")
                    if(column %in% c("scoreSR")){
                      Gviz::displayPars(dt)$ylim <- c(0,1)
                    } else if(!is.null(args[["ylim"]])){
                      Gviz::displayPars(dt)$ylim <- args[["ylim"]]
                    }
                    Gviz::displayPars(dt)$background.title <- "#FFFFFF"
                    Gviz::displayPars(dt)$fontcolor.title <- "#000000"
                    Gviz::displayPars(dt)$col.axis <- "#000000"
                    if(!is.null(args[["data.track.pars"]])){
                      Gviz::displayPars(dt) <- args[["data.track.pars"]]
                    }
                    dt
                  })
    names(dts) <- colnames(mcols(data))
    dts
  }
)

#' @rdname ModAlkAnilineSeq
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModSetAlkAnilineSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("scoreNC","scoreSR","ends"),
                        window.size = 15L, ...) {
    type <- match.arg(type,c("scoreNC","scoreSR","ends"))
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModAlkAnilineSeq
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModSetAlkAnilineSeq"),
  definition = function(x, name, from, to,
                        type = c("scoreNC","scoreSR","ends"),
                        ...) {
    type <- match.arg(type,c("scoreNC","scoreSR","ends"))
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)
