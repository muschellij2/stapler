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
#' @param set_origin Should the origin be set to the same if x is a
#' set of images, including \code{niftiImage}s.
#' @export
staple = function(
  x,
  ...,
  set_origin = FALSE) {
  UseMethod("staple")
}

#' @rdname staple
#' @export
staple.default = function(
  x,
  ...,
  set_origin = FALSE){
  res = staple_multi_mat(x, ...)
  return(res)
}


#' @rdname staple
#' @export
staple.list = function(
  x,
  ...,
  set_origin = FALSE){

  res = staple_multi_img(x, set_origin = set_origin, ...)
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
  set_origin = FALSE){

  res = staple_multi_img(x, set_origin = set_origin, ...)
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
staple.array = function(
  x,
  ...,
  set_origin = FALSE){
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


