#' Run STAPLE on a set of nifti images
#'
#' @param x Character vector of filenames or list of arrays/images
#' @param set_orient Should the orientation be set to the same if the images are
#' \code{niftiImage}s
#' @param verbose print diagnostic messages
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
#' x = lapply(seq(n), function(i) {
#'    x = rbinom(n = r, size = 1, prob = 0.5)
#'    array(x, dim = c(10,10, 10))
#'  })
#' staple_out = staple_bin_img(x, set_orient = FALSE)
#' res = staple(x)
#' testthat::expect_equal(staple_out$sensitivity,
#' res$sensitivity)
#'
#' @importFrom RNifti readNifti niftiHeader updateNifti orientation
#' @importFrom RNifti "orientation<-"
staple_bin_img = function(
  x,
  set_orient = FALSE,
  verbose = TRUE,
  ...) {

  if (verbose) {
    message("Reshaping images")
  }
  x = reshape_img(x = x, set_orient = set_orient,
                  verbose = verbose)
  first_image = x$first_image
  all_nifti = x$all_nifti
  x = x$x
  if (verbose) {
    message("Running STAPLE for binary matrix")
  }
  res = staple_bin_mat(x, verbose = verbose, ...)

  if (verbose) {
    message("Creating output probability image/array")
  }
  outimg = array(res$probability,
                 dim = dim(first_image))
  if (all_nifti) {
    hdr = RNifti::niftiHeader(first_image)
    hdr$cal_max = 1
    hdr$cal_min = 0
    hdr$datatype = 16
    hdr$bitpix = 32

    outimg = RNifti::updateNifti(
      outimg, template = hdr)
  }

  if (verbose) {
    message("Creating output prior image/array")
  }
  priorimg = array(res$prior,
                   dim = dim(first_image))
  if (all_nifti) {
    priorimg = RNifti::updateNifti(
      priorimg, template = hdr)
  }

  if (verbose) {
    message("Creating label image (probability >= 0.5)")
  }
  label = array(res$label,
                dim = dim(first_image))
  if (all_nifti) {
    hdr$datatype = 2
    hdr$bitpix = 8

    label = RNifti::updateNifti(
      label, template = hdr)
  }


  res$probability = outimg
  res$label = label
  res$prior = priorimg
  rm(list = c("outimg", "priorimg", "label"))
  gc()
  return(res)

}

#' @export
#' @rdname staple_bin_img
#' @examples
#' n = 5
#' r = 1000
#' x = lapply(seq(n), function(i) {
#'    x = rbinom(n = r, size = 5, prob = 0.5)
#'    array(x, dim = c(10,10, 10))
#'  })
#' staple_out = staple_multi_img(x, set_orient = FALSE)
staple_multi_img = function(
  x,
  set_orient = FALSE,
  verbose = TRUE,
  ...) {

  if (verbose) {
    message("Reshaping images")
  }
  x = reshape_img(x = x, set_orient = set_orient,
                  verbose = verbose)
  first_image = x$first_image
  all_nifti = x$all_nifti
  x = x$x
  res = staple_multi_mat(x, ...)

  if (all_nifti) {
    hdr = RNifti::niftiHeader(first_image)
    hdr$cal_max = 1
    hdr$cal_min = 0
    hdr$datatype = 16
    hdr$bitpix = 32
  }

  if (verbose) {
    message("Creating output probability images/arrays")
  }
  n_level = ncol(res$probability)
  outimg = lapply(seq(n_level), function(ind) {
    probability = res$probability[, ind]
    outimg = array(
      probability,
      dim = dim(first_image))
    if (all_nifti) {
      outimg = RNifti::updateNifti(
        outimg, template = hdr)
    }
    return(outimg)
  })
  names(outimg) = colnames(res$probability)
  res$probability = outimg
  rm(list = "outimg"); gc()

  priorimg = lapply(seq(n_level), function(ind) {
    probability = res$prior[, ind]
    outimg = array(
      probability,
      dim = dim(first_image))
    if (all_nifti) {
      outimg = RNifti::updateNifti(
        outimg, template = hdr)
    }
    return(outimg)
  })
  names(priorimg) = colnames(res$prior)
  res$prior = priorimg
  rm(list = "priorimg"); gc()

  if (verbose) {
    message("Creating output label image/array")
  }
  label = array(
    res$label,
    dim = dim(first_image))

  if (all_nifti) {
    hdr$datatype = 8
    hdr$bitpix = 32
    label = RNifti::updateNifti(
      label, template = hdr)
  }
  res$label = label

  return(res)
}

