#' @include RNAmodR.AlkAnilineSeq.R
#' @include Modifier-AlkAnilineSeq-class.R
NULL

RNAMODR_AAS_PLOT_DATA <- c("ends",
                           "scoreNC",
                           "scoreSR")
RNAMODR_AAS_PLOT_DATA_DEFAULT <- c("scoreNC")
RNAMODR_AAS_PLOT_DATA_NAMES <- c(ends = "5'-ends",
                                 scoreNC = "Normalized Cleavage",
                                 scoreSR = "Stop ratio")
RNAMODR_AAS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   scoreNC = "#DECBE4",
                                   scoreSR = "#CCEBC5")

.norm_viz_mod_aas_args <- function(input, type){
  if(!all(type %in% RNAMODR_AAS_PLOT_DATA)){
    stop("Type '",type,"' is not valid. Valid types are: '",
         paste0(RNAMODR_AAS_PLOT_DATA, collapse = "','"),"'.",
         call. = FALSE)
  }
  colour <- input[["colour"]]
  if(!is.null(input[["colour"]])){
    colour <- RNAmodR:::.norm_viz_colour(input[["colour"]], type)
  } else {
    colour <- RNAMODR_AAS_PLOT_DATA_COLOURS[type]
  }
  input <- list(type = type,
                colour = colour)
  input
}

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = function(x, name, type, ...) {
    args <- .norm_viz_mod_aas_args(list(...), type)
    data <- RNAmodR:::.get_data_for_visualization(x, name)
    data <- unlist(data)
    dts <- lapply(
      args[["type"]],
      function(t){
        dt <- Gviz::DataTrack(range = data[,t],
                              groups = t,
                              name = RNAMODR_AAS_PLOT_DATA_NAMES[t],
                              col = args[["colour"]][t],
                              type = "histogram")
        Gviz::displayPars(dt)$background.title <- "#FFFFFF"
        Gviz::displayPars(dt)$fontcolor.title <- "#000000"
        Gviz::displayPars(dt)$col.axis <- "#000000"
        Gviz::displayPars(dt) <- args[names(args) != "type"]
        dt
      })
    names(dts) <- args[["type"]]
    dts
  }
)

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModAlkAnilineSeq",
                        coord = "GRanges"),
  definition = function(x, coord, type = c("ends","scoreNC","scoreSR"),
                        window.size = 15L, ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,RNAMODR_AAS_PLOT_DATA, several.ok = TRUE)
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModAlkAnilineSeq"),
  definition = function(x, name, from, to, type = c("ends","scoreNC","scoreSR"),
                        ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,RNAMODR_AAS_PLOT_DATA, several.ok = TRUE)
    callNextMethod(x = x, name = name, from = from, to = to, type = type, ...)
  }
)

#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModSetAlkAnilineSeq",
                        coord = "GRanges"),
  definition = function(x, coord, type = c("scoreNC","scoreSR","ends"),
                        window.size = 15L, ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,RNAMODR_AAS_PLOT_DATA, several.ok = TRUE)
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModAlkAnilineSeq-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModSetAlkAnilineSeq"),
  definition = function(x, name, from, to, type = c("scoreNC","scoreSR","ends"),
                        ...) {
    if(missing(type)){
      type <- RNAMODR_AAS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,RNAMODR_AAS_PLOT_DATA, several.ok = TRUE)
    callNextMethod(x = x, name = name, from = from, to = to, type = type, ...)
  }
)
