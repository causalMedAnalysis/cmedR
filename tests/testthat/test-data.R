test_that("packaged data sets load and have expected structure", {
  data("nlsy", package = "cmedR", envir = environment())
  data("brader_clean", package = "cmedR", envir = environment())

  expect_s3_class(nlsy, "data.frame")
  expect_gt(nrow(nlsy), 0)
  expect_true(all(
    c("cesd_age40", "ever_unemp_age3539", "log_faminc_adj_age3539", "att22",
      "female", "black", "hispan", "paredu", "parprof", "parinc_prank",
      "famsize", "afqt3") %in% names(nlsy)
  ))

  expect_s3_class(brader_clean, "data.frame")
  expect_gt(nrow(brader_clean), 0)
  expect_true(all(
    c("immigr", "tone_eth", "p_harm", "emo") %in% names(brader_clean)
  ))
})
