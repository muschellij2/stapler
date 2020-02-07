#' Generic STAPLE Algorithm
#'
#' Tries to do the correct STAPLE algorithm (binary/multi-class) for the
#' type of input (array/matrix/list of images/filenames of images)
#'
#' @param x a nxr matrix where there are n raters and r elements rated,
#' a list of images, or a character vector.  Note, \code{\link{readNifti}}
#' is used for image filenames
#' @param ... Options for STAPLE, see \code{\link{staple_bin_mat}} and
#' \code{\link{staple_multi_mat}}
#' @param set_orient Should the orientation be set to the same if x is a
#' set of images, including \code{niftiImage}s.
#' @export
#' @examples
#' n = 5
#' r = 1000
#' sens = c(0.8, 0.9, 0.8, 0.5, 0.8)
#' spec = c(0.9, 0.75, 0.99, 0.98, 0.92)
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(20171120)
#' n_1 = 200
#' n_0 = r - n_1
#' truth = c(rep(0, n_0), rep(1, n_1))
#' pred_1 = rbinom(n = n, size = n_1, prob = sens)
#' pred_0 = rbinom(n = n, size = n_0, prob = spec)
#' pred_0 = sapply(pred_0, function(n) {
#'    sample(c(rep(0, n), rep(1, n_0 -n)))
#' })
#' pred_1 = sapply(pred_1, function(n) {
#'    sample(c(rep(1, n), rep(0, n_1 -n)))
#' })
#' pred = rbind(pred_0, pred_1)
#' true_sens = colMeans(pred[ truth == 1, ])
#' true_spec = colMeans(1-pred[ truth == 0, ])
#' x = t(pred)
#' staple_out = staple(x)
#' print(staple_out$sensitivity)
#' testthat::expect_equal(staple_out$sensitivity[,"1"],
#' c(0.781593858553476, 0.895868301462594,
#' 0.760514086161722, 0.464483444340873,
#' 0.765239314719065))
#' staple_out_prior = staple(x, prior = rep(0.5, r))
#'
#' testthat::expect_equal(staple_out_prior$sensitivity[, "1"],
#' c(0.683572080864211, 0.821556768891859,
#' 0.619166852992802, 0.389409921992467, 0.67042085955546))
#'
#' res_bin = staple_bin_mat(x, prior = rep(0.5, 1000))
#' testthat::expect_equal(staple_out_prior$sensitivity[,"1"],
#' res_bin$sensitivity)
staple = function(
  x,
  ...,
  set_orient = FALSE) {
  UseMethod("staple")
}

#' @rdname staple
#' @export
staple.default = function(
  x,
  ...,
  set_orient = FALSE){
  res = staple_multi_mat(x, ...)
  return(res)
}


#' @rdname staple
#' @export
staple.list = function(
  x,
  ...,
  set_orient = FALSE){

  res = staple_multi_img(x, set_orient = set_orient, ...)
  if (length(res$probability) == 2) {
    res$probability = res$probability[[2]]
    res$prior = res$prior[[2]]
    res$sensitivity = res$sensitivity[, 2]
    res$specificity = res$specificity[, 2]
  }
  return(res)
}

#' @rdname staple
#' @export
staple.character = function(
  x,
  ...,
  set_orient = FALSE){

  res = staple_multi_img(x, set_orient = set_orient, ...)
  if (length(res$probability) == 2) {
    res$probability = res$probability[[2]]
    res$prior = res$prior[[2]]
    res$sensitivity = res$sensitivity[, 2]
    res$specificity = res$specificity[, 2]
  }
  return(res)
}

# array covers array/nifti/matrix
#' @rdname staple
#' @export
#' @examples
#' n = 5
#' r = 1000
#' x = lapply(seq(n), function(i) {
#'    x = rbinom(n = r, size = 1, prob = 0.5)
#'    array(x, dim = c(10,10, 10))
#'  })
#' mat = sapply(x, c)
#' staple_out = staple_bin_img(x, set_orient = FALSE)
#' res_mat = staple(t(mat))
#' testthat::expect_equal(staple_out$sensitivity, res_mat$sensitivity[, "1"])
staple.array = function(
  x,
  ...,
  set_orient = FALSE){
  ndim = length(dim(x))
  if (ndim > 2) {
    stop(paste0("Array of more than 2 dimensions given.  ",
                "If 4D array (such as images), convert to list of ",
                "3D images"))
  }

  res = staple_multi_mat(x, ...)
  if (ncol(res$probability) == 2) {
    res$probability = res$probability[,2]
    res$prior = res$prior[,2]
    res$sensitivity = res$sensitivity[, 2]
    res$specificity = res$specificity[, 2]
  }
  return(res)
}


