testthat::test_that("Staple binary matrix", {
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(20171120)
  n = 5
  r = 1000
  sens = c(0.8, 0.9, 0.8, 0.5, 0.8)
  spec = c(0.9, 0.75, 0.99, 0.98, 0.92)
  n_1 = 200
  n_0 = r - n_1
  truth = c(rep(0, n_0), rep(1, n_1))
  pred_1 = rbinom(n = n, size = n_1, prob = sens)
  pred_0 = rbinom(n = n, size = n_0, prob = spec)
  pred_0 = sapply(pred_0, function(n) {
    sample(c(rep(0, n), rep(1, n_0 - n)))
  })
  pred_1 = sapply(pred_1, function(n) {
    sample(c(rep(1, n), rep(0, n_1 - n)))
  })
  pred = rbind(pred_0, pred_1)
  true_sens = colMeans(pred[ truth == 1, ])
  true_spec = colMeans(1 - pred[ truth == 0, ])
  x = t(pred)

  # need test for getRversion() >= numeric_version("3.6.0")
  testthat::expect_message({res = staple_bin_mat(x)})
  testthat::expect_equal(
    res$sensitivity,
    c(0.781593858553476, 0.895868301462594,
      0.760514086161722, 0.464483444340873,
      0.765239314719065))
  testthat::expect_equal(
    res$specificity,
    c(0.902896626562703, 0.770915583628547, 0.994826018925032,
      0.979916045260754,
      0.935564390547129))
  table(res$label, truth)
  accuracy = mean(res$label == truth)
  testthat::expect_equal(accuracy, 0.974)

  testthat::expect_silent({
    res2 = staple_bin_mat(x, prior = rep(0.5, r),
                          verbose = FALSE)
  })
  testthat::expect_equal(res2$sensitivity,
                         c(0.683572080864211, 0.821556768891859,
                           0.619166852992802, 0.389409921992467,
                           0.67042085955546)
  )
  testthat::expect_equal(res2$specificity,
                         c(0.919431705156552, 0.795365633127613,
                           0.99999999994706, 0.986146453778281,
                           0.95467606186872)
  )
  table(res2$label, truth)

  #######################################
  # Given only 2 classes - should give same
  #######################################
  testthat::expect_message({
    multi_res = staple_multi_mat(x)
  })
  testthat::expect_equal(res$label*1, multi_res$label*1)
  testthat::expect_equal(res$sensitivity, multi_res$sensitivity[, "1"])
  testthat::expect_equal(res$specificity, multi_res$specificity[, "1"])



})
