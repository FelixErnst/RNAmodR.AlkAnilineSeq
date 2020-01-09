#' @title RNAmodR.AlkAnilineSeq
#' 
#' @author Felix G M Ernst [aut], Denis L J Lafontaine [fnd]
#' 
#' @description
#' `RNAmodR.AlkAnilineSeq` implements the detection of 7-methyl guanosine,
#' 3-methyl cytidine and dihydrouridine from AlkAnilineseq data using the 
#' workflow and class the package `RNAmodR` provides.
#' 
#' @seealso Further details are described in the man pages of the 
#' \code{\link[RNAmodR:Modifier-class]{Modifier}} object and the vignettes.
#'
#' @docType package
#' @name RNAmodR.AlkAnilineSeq
NULL

#' @import methods
#' @import RNAmodR
#' @import S4Vectors
#' @import BiocGenerics
#' @import IRanges
#' @import GenomicRanges
#' @import Gviz
NULL

#' @name RNAmodR.AlkAnilineSeq-datasets
#' @title Example data in the RNAmodR.AlkAnilineSeq package
#' @description 
#' This contains an example ModifierSet object of type ModSetAlkAnilineSeq
#' @docType data
#' @format a \code{ModSetAlkAnilineSeq} instance
#' @keywords datasets
#' 
#' @usage data(msaas)
"msaas"