context("AlkAnilineSeq")
# .get_heatmap_data
test_that("AlkAnilineSeq:",{
  # argument normalization
  input <- list()
  actual <- RNAmodR.AlkAnilineSeq:::.norm_aas_args(input)
  expect_type(actual,"list")
  expect_named(actual,c("minCoverage",
                        "minReplicate",
                        "findMod",
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
})
