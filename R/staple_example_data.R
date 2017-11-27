#' STAPLE Example Data
#'
#' @return Character vector of filenames
#' @export
#'
#' @examples
#' staple_example_data()
staple_example_data = function() {
  fnames = paste0("seg_", 1:3, ".nii.gz")
  system.file("extdata", fnames, package = "stapler")
}
