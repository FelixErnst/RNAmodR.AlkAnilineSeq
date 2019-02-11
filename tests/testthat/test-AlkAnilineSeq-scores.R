context("AlkAnilineSeq base score")
# .get_heatmap_data
test_that("AlkAnilineSeq base score",{
  # base score selection
  df <- S4Vectors::DataFrame(means.treated.G = c("100","0"),
                             means.treated.A = c("0","100"),
                             means.treated.T = c("0","100"),
                             means.treated.C = c("0","100"))
  x <- IRanges::SplitDataFrameList(df,df)
  seq <- Biostrings::RNAStringSet(c("AG","CU"))
  actual <- RNAmodR.AlkAnilineSeq:::.calculate_base_score(x,seq)
  expect_equal(actual,c(0,0,0,100))
})
