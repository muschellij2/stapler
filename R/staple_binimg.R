ensure_Nifti = function(x) {
  if (is.list(x)) {
    x = lapply(x, ensure_Nifti)
    return(x)
  }
  if (is.character(x)) {
    if (length(x) > 0) {
      x = lapply(x, readNifti)
    } else {
      x = readNifti(x)
    }
  }
  # allow for arrays
  if (is.array(x)) {
    return(x)
  }
  if (!inherits(x, "niftiImage")) {
    warning("x is not a RNifti::niftiImage image - may not work!")
  }
  x
}

reshape_img = function(x, set_origin = TRUE) {
  x = ensure_Nifti(x)
  # orientations = sapply(imgs, orientation)
  first_image = x[[1]]

  if (set_origin) {
    ori = orientation(first_image)
    x = lapply(x,
               function(x) {
                 orientation(x) = ori
                 x
               })
  }

  x = t(sapply(x, c))
  L = list(x = x,
           first_image = first_image)
  return(L)
}

#' Run STAPLE on a set of nifti images
#'
#' @param x Character vector of filenames or list of arrays/images
#' @param set_origin Should the origin be set to the same if they images are
#' \code{niftiImage}s
#' @param ... Additional arguments to \code{\link{staple_bin_mat}}
#'
#' @return A list similar to \code{\link{staple_bin_mat}}, but
#' has a resulting image
#' @export
#'
#' @examples
#'
#' @examples
#' n = 5
#' r = 1000
#' sens = c(0.8, 0.9, 0.8, 0.5, 0.8)
#' spec = c(0.9, 0.75, 0.99, 0.98, 0.92)
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
#' x = lapply(seq(n), function(x) {
#'    array(pred[, x], dim = c(10,10, 10))
#'  })
#' staple_out = staple_bin_img(x, set_origin = FALSE)
#'
#' @importFrom RNifti readNifti dumpNifti updateNifti orientation
#' @importFrom RNifti "orientation<-"
staple_bin_img = function(
  x,
  set_origin = FALSE,
  ...) {

  x = reshape_img(x = x, set_origin = set_origin)
  first_image = x$first_image
  x = x$x
  res = staple_bin_mat(x, ...)


  outimg = array(res$probability,
                 dim = dim(first_image))

  hdr = RNifti::dumpNifti(first_image)
  hdr$cal_max = 1
  hdr$cal_min = 0
  hdr$datatype = 16

  outimg = RNifti::updateNifti(
    outimg, template = hdr)

  res$outimg = outimg
  return(res)

}

#' @export
#' @rdname staple_bin_img
staple_multi_img = function(
  x,
  set_origin = FALSE,
  ...) {

  x = reshape_img(x = x, set_origin = set_origin)
  first_image = x$first_image
  x = x$x
  res = staple_multi_mat(x, ...)

  prob = res$probability
  n_level = ncol(prob)
  outimg = lapply(seq(n_level), function(ind) {
    probability = prob[, ind]
    outimg = array(
      probability,
      dim = dim(first_image))

    hdr = RNifti::dumpNifti(first_image)
    hdr$cal_max = 1
    hdr$cal_min = 0
    hdr$datatype = 16

    outimg = RNifti::updateNifti(
      outimg, template = hdr)
  })
  res$outimg = outimg
  return(res)
}
