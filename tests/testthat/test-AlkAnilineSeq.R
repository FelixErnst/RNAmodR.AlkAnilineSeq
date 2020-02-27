context("AlkAnilineSeq")
# .get_heatmap_data
test_that("AlkAnilineSeq:",{
  # argument normalization
  input <- list()
  actual <- RNAmodR.AlkAnilineSeq:::.norm_aas_args(input)
  expect_type(actual,"list")
  expect_named(actual,c("minCoverage",
                        "minReplicate",
                        "find.mod",
                        "minLength",
                        "minSignal",
                        "minScoreNC",
                        "minScoreSR",
                        "minScoreBaseScore",
                        "scoreOperator"))
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minCoverage = 10)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minCoverage = 10L)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minReplicate = 1)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minReplicate = 1L)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minSignal = 10)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minSignal = 10L)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minLength = 9)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minLength = 9L)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minSignal = 10)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minSignal = 10L)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minScoreNC = "a")))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minScoreNC = 50)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minScoreSR = 2L)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(minScoreSR = 0.5)),
               actual)
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(scoreOperator = 2L)))
  expect_equal(RNAmodR.AlkAnilineSeq:::.norm_aas_args(list(scoreOperator = "&")),
               actual)
  data <- data.frame(scoreNC = "1", scoreSR = "a",
                     stringsAsFactors = FALSE)
  actual <- RNAmodR.AlkAnilineSeq:::.get_aas_scores(data)
  expect_type(actual,"list")
  expect_named(actual,c("score","scoreSR"))
  # settings
  data("msaas", package = "RNAmodR.AlkAnilineSeq")
  expect_type(settings(msaas),"list")
  expect_type(settings(msaas[[1]]),"list")
  expect_error(settings(msaas[[1]]) <- 1,
               "'value' has to be a named.")
  expect_error(settings(msaas[[1]]) <- c(1),
               "'value' has to be a named.")
  # settings(msaas[[1]]) <- c(minCoverage = 11L)
  # expect_equal(settings(msaas[[1]])$minCoverage, 11L)
  # apply_Np1_offset 
  mod <- aggregate(sequenceData(msaas[[1]]), condition = "Treated")$NormEnd5SequenceData
  actual <- RNAmodR.AlkAnilineSeq:::.apply_Np1_offset(mod)
  expect_named(actual)
  expect_equal(mod[[1]][5,1], 13L)
  expect_equal(actual[[1]][5,1], 20L)
  # .aggregate_aas 
  actual <- RNAmodR.AlkAnilineSeq:::.aggregate_aas(msaas[[1]]) 
  expect_named(actual)
  expect_s4_class(actual,"CompressedSplitDFrameList")
  expect_equivalent(actual[[1]][1,3], 0.1046218)
  expect_equal(colnames(actual[[1]]),c("ends","scoreNC","scoreSR","baseScore"))
  # findMod
  actual <- findMod(msaas[[1]])
  expect_s4_class(actual,"GRanges")
  expect_equal(colnames(mcols(actual)),c("mod","source","type","score",
                                         "scoreSR","Parent"))
  expect_equal(unique(mcols(actual)$mod),c("m7G","D","m3C"))
  expect_equal(length(actual),9L)
  # show
  # expect_output(expect_warning(show(msaas),
                               # "Settings were changed after data aggregation or"))
  # expect_warning(modifications(msaas),
                 # "Settings were changed after data aggregation or")
  expect_true(all(modifications(msaas)[[1]] == actual))
  # plotting
  expect_error(RNAmodR.AlkAnilineSeq:::.norm_viz_mod_aas_args(list(), "1"),
               "Type '1' is not valid")
  actual <- RNAmodR.AlkAnilineSeq:::.norm_viz_mod_aas_args(list(), "scoreNC")
  expect_type(actual,"list")
  expect_named(actual,c("type","colour"))
  expect_equal(actual[["type"]],"scoreNC")
  actual <- RNAmodR.AlkAnilineSeq:::.norm_viz_mod_aas_args(list(), "scoreSR")
  expect_equal(actual[["type"]],"scoreSR")
  actual <- RNAmodR.AlkAnilineSeq:::.norm_viz_mod_aas_args(list(), "ends")
  expect_equal(actual[["type"]],"ends")
  # getDataTrack
  expect_error(getDataTrack(msaas[[1]]),
               'argument "type" is missing, with no default')
  expect_error(getDataTrack(msaas[[1]],type="score"),
               "Type 'score' is not valid")
  expect_error(getDataTrack(msaas[[1]],type="scoreNC"),
               'argument "name" is missing')
  actual <- getDataTrack(msaas[[1]],name="1",type="scoreNC")
  expect_type(actual,"list")
  expect_named(actual,"scoreNC")
  expect_s4_class(actual[[1]],"DataTrack")
  # plotData
  actual <- plotData(msaas[[1]],name="1",type="scoreNC")
  expect_type(actual,"list")
  expect_named(actual,c("Normalized Cleavage","SequenceTrack","titles"))
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[2]],"SequenceRNAStringSetTrack")
  expect_s4_class(actual[[3]],"ImageMap")
  actual <- plotData(msaas,name="1",type="scoreNC")
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[2]],"DataTrack")
  expect_s4_class(actual[[3]],"DataTrack")
  expect_s4_class(actual[[4]],"SequenceRNAStringSetTrack")
  expect_s4_class(actual[[5]],"ImageMap")
})
