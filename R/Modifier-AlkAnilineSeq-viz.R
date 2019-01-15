#' @include RNAmodR.AlkAnilineSeq.R
#' @include Modifier-AlkAnilineSeq-class.R
NULL

RNAMODR_AAS_PLOT_DATA <- c("ends",
                           "score")
RNAMODR_AAS_PLOT_DATA_NAMES <- c(ends = "5'-ends",
                                 score = "Score AlkAnilineSeq")
RNAMODR_AAS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   score = "#DECBE4")
#' @rdname ModAlkAnilineSeq
#' 
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
  definition = function(x,
                        coord,
                        type = c("ends","score"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA
    }
    callNextMethod(x = x,
                   coord = coord,
                   type = type,
                   window.size = window.size,
                   ...)
  }
)
#' @rdname ModAlkAnilineSeq
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = function(x,
                        name,
                        from,
                        to,
                        type = c("ends","score"),
                        ...) {
    browser()
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA
    }
    callNextMethod(x = x,
                   name,
                   from,
                   to,
                   type = type,
                   ...)
  }
)

setMethod(
  f = ".dataTracks",
  signature = signature(x = "ModAlkAnilineSeq",
                        data = "GRanges",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x,
                        data,
                        seqdata,
                        sequence,
                        args) {
    requireNamespace("Gviz")
    browser()
    n <- ncol(mcols(data))
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_AAS_PLOT_DATA_COLOURS
    }
    dts <- lapply(seq_len(n),
                  function(i){
                    column <- colnames(mcols(data)[i])
                    colour <- colour[column]
                    name <- RNAMODR_AAS_PLOT_DATA_NAMES[column]
                    dt <- DataTrack(data,
                              data = column,
                              name = name,
                              fill = colour,
                              type = "histogram")
                    if(column %in% c("score")){
                      displayPars(dt)$ylim = c(0,1)
                    }
                    displayPars(dt)$background.title <- "#FFFFFF"
                    displayPars(dt)$fontcolor.title <- "#000000"
                    displayPars(dt)$col.axis <- "#000000"
                    dt
                  })
    names(dts) <- colnames(mcols(data))
    dts
  }
)
