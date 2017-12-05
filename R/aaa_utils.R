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
    return(x)
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

reshape_img = function(x, set_orient = TRUE, verbose = TRUE) {
  x = ensure_Nifti(x)
  msg = "Using a staple_*_img function, but images not passed!"
  if (is.matrix(x)) {
    stop(msg)
  }
  if (!is.list(x)) {
    stop(msg)
  }
  # orientations = sapply(imgs, orientation)
  first_image = x[[1]]

  nifti_images = sapply(x, inherits, what = "niftiImage")
  all_nifti = all(nifti_images)
  any_nifti = any(nifti_images)
  if (any_nifti & !all_nifti) {
    stop("Some inputs are niftiImages and others are not, failing")
  }
  if (all_nifti && verbose) {
    message("All images are niftiImage ojects")
  }
  if (set_orient) {
    if (all_nifti) {
      ori = orientation(first_image)
      x = lapply(x,
                 function(x) {
                   orientation(x) = ori
                   x
                 })
    } else {
      warning("set_orient = TRUE, but niftiImages not passed!")
    }
  }
  dims = lapply(x, dim)
  res = sapply(dims,
               identical,
               x = dims[[1]])
  res = all(res)
  if (!all(res)) {
    stop("Images are not all the same dimensions!")
  }
  x = t(sapply(x, c))
  L = list(x = x,
           first_image = first_image,
           all_nifti = all_nifti)
  return(L)
}
