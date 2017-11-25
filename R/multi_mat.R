#' STAPLE on Multi-class matrix
#'
#' @param x a nxr matrix where there are n raters and r elements rated
#' @param sens_init Initalize matrix for sensitivity (p)
#' @param spec_init  Initalize matrix for specificity (q)
#' @param max_iter Maximum number of iterations to run
#' @param tol Tolerance for convergence
#' @param prior Either "mean" or a matrix of prior probabilities,
#' @param verbose print diagnostic messages
#' @param trace Number for modulus to print out verbose iterations
#' @param ties.method Method passed to \code{\link{max.col}}
#' for hard segmentation
#'
#' @return List of matrix output sensitivities, specificities, and
#' matrix of probabilities
#' @export
#'
#' @examples
#' rm(list = ls())
#' x = matrix(rbinom(5000, size = 5, prob = 0.5), ncol = 1000)
#'   sens_init = 0.99999
#'   spec_init = 0.99999
#'   max_iter = 10000
#'   tol = .Machine$double.eps
#'   prior = "mean"
#'   verbose = TRUE
#'   trace = 25
#'   ties.method = "first"
#'
#' res = staple_multi_mat(x)
#'
#' @importFrom matrixStats colProds colVars
staple_multi_mat = function(
  x,
  sens_init = 0.99999,
  spec_init = 0.99999,
  max_iter = 10000,
  tol = .Machine$double.eps,
  prior = "mean",
  verbose = TRUE,
  trace = 25,
  ties.method = c("first", "random", "last")
){

  n_readers = nrow(x)
  n_all_voxels = ncol(x)

  if (n_readers > n_all_voxels) {
    warning(paste0(
      "Number of readers larger than number of elements.",
      "Are you sure matrix x is nxr?")
    )
  }

  stopifnot(!any(is.na(x)))
  umat = sort(unique(c(x)))
  umat = as.numeric(umat)
  n_levels = length(umat)

  if (verbose) {
    message("Making multiple, matrices. Hot-one encode")
  }
  xmats = lapply(umat, function(val) {
    x == val
  })
  names(xmats) = umat

  # not_all_same = lapply(xmats, function(x) {
  #   cs = colSums(x)
  #   all_zero = cs == 0
  #   # only_one = cs == 1
  #   # if all vote yes - then yes
  #   all_one = cs == n_readers
  #   not_all_same = !all_zero & !all_one
  #
  #   stopifnot(!any(is.na(not_all_same)))
  #   return(list(not_all_same = not_all_same, all_one = all_one))
  # })
  # not_all_same = Reduce("|", not_all_same)
  not_all_same = matrixStats::colVars(x) > 0

  ####################################
  # Keeping only voxels with more than 1 says yes
  ####################################
  # prior = match.arg(prior)
  p1 = prior[1]
  if (p1 == "mean") {
    f_t_i = sapply(xmats, colMeans, na.rm = TRUE)
    prior = f_t_i
  } else {
    stop("not implemented")
    prior = as.matrix(prior)
    n_prior = ncol(prior)
    if (n_prior != n_all_voxels) {
      stop("Prior does not have same number of rated elements!")
    }
    stopifnot(!any(is.na(prior)))
    f_t_i = prior
    ####################
    # wrong
    if (any(prior %in% c(0, 1))) {
      warning("Some elements in prior are in {0, 1}")
    }
    # all_one = all_one | prior == 1
    # all_zero = all_zero | prior == 0
    # not_all_same = !all_zero & !all_one
    # stopifnot(!any(is.na(not_all_same)))
  }
  if (verbose) {
    message("Removing elements where all raters agree")
  }
  # mats = lapply(xmats, function(r) r[, not_all_same])
  # mats = xmats
  mats = lapply(xmats, function(r) r[, not_all_same])
  f_t_i = f_t_i[not_all_same,]

  d_f_t_i = 1 - f_t_i


  # n_voxels = sum(not_all_same)
  n_voxels = ncol(mats[[1]])
  dmats = lapply(
    mats, function(mat) {
      # dmat = (1L - mat) > 0
      # class(dmat) = "logical"
      # dmat[dmat == 0] = NA
      dmat = !mat
      return(dmat)
    })

  # doing this for na.rm arguments
  # mats = lapply(
  #   mats, function(mat) {
  #     mat[ mat == 0] = NA
  #     mat
  #   })

  ###################
  #initialize
  p = matrix(sens_init, nrow = n_readers, ncol = n_levels)
  q = matrix(spec_init, nrow = n_readers, ncol = n_levels)


  eps = sqrt(tol)

  # mat is D
  ### run E Step
  for (iiter in seq(max_iter)) {
    # pmat = sapply(seq(n_levels), function(ind) {
    #   p = p[, ind]
    #   mat = mats[[ind]]
    #   pmat = p * mat
    #   pmat = matrixStats::colProds(pmat, na.rm = TRUE)
    # })
    # sep_pmat = sapply(seq(n_levels), function(ind) {
    #   p = p[, ind]
    #   dmat = dmats[[ind]]
    #   sep_pmat = (1 - p) * dmat
    #   sep_pmat =  matrixStats::colProds(sep_pmat, na.rm = TRUE)
    # })
    #
    # qmat = sapply(seq(n_levels), function(ind) {
    #   q = q[, ind]
    #   mat = mats[[ind]]
    #   qmat = q * mat
    #   qmat =  matrixStats::colProds(qmat, na.rm = TRUE)
    # })
    # sep_qmat = sapply(seq(n_levels), function(ind) {
    #   q = q[, ind]
    #   dmat = dmats[[ind]]
    #   sep_qmat = (1 - q) * dmat
    #   sep_qmat =  matrixStats::colProds(sep_qmat, na.rm = TRUE)
    # })

    W_i = sapply(seq(n_levels), function(ind) {
      p = p[, ind]
      q = q[, ind]
      ft = f_t_i[, ind]
      dft = d_f_t_i[, ind]

      mat = mats[[ind]]
      dmat = dmats[[ind]]

      # E Step
      a_i = p ^ mat * (1 - p) ^ dmat
      a_i = ft * matrixStats::colProds(a_i)

      b_i = q ^ dmat * (1 - q) ^ mat
      b_i = dft * matrixStats::colProds(b_i)
      W_i = a_i / (a_i + b_i)
      W_i
    })

    W_i = W_i / rowSums(W_i)

    W_is = lapply(seq(n_levels), function(ind) {
      W_i[, ind]
    })

    ##########################
    # do these make sense to do?
    ##########################
    # W_i = pmin(W_i, 1 - eps)
    # W_i = pmax(W_i, eps)
    f = mapply(function(mat, dmat, W_i) {
      sum_w = sum(W_i)

      new_p  = t(mat) * W_i
      # new_p  = colSums(new_p,	na.rm = TRUE)
      new_p  = colSums(new_p)
      new_p = new_p/(sum_w + eps)

      new_q  = t(dmat) * (1 - W_i)
      # new_q  = colSums(new_q,	na.rm = TRUE)
      new_q  = colSums(new_q)
      new_q = new_q/(n_voxels - sum_w + eps)
      return(list(new_p = new_p, new_q = new_q))
    }, mats, dmats, W_is, SIMPLIFY = FALSE)
    new_p = sapply(f, function(x){
      x$new_p
    })
    new_q = sapply(f, function(x){
      x$new_q
    })

    diff_p = abs(p - new_p)
    diff_q = abs(q - new_q)
    diff = max(c(diff_p, diff_q))
    if (diff <= tol) {
      if (verbose) {
        message("Convergence!")
      }
      break
    } else {
      if (verbose) {
        if (iiter %% trace == 0) {
          message(paste0("iter: ", iiter,
                         ", diff: ", diff))
        }
      }
    }

    p = new_p
    q = new_q
  }

  colnames(W_i) = umat
  stopifnot(!any(is.na(W_i)))

  outimg = matrix(NA, nrow = n_all_voxels, ncol = n_levels)
  colnames(outimg) = umat
  outimg[not_all_same, ] = W_i
  # those all the same are given prob 1 of that class
  xind = x[1, !not_all_same]
  sub = outimg[!not_all_same, ]
  sub_replacement = sapply(umat, function(r) r == xind)
  stopifnot(all(dim(sub) == dim(sub_replacement)))
  # sub[, umat == xind] = 1
  # sub[, umat != xind] = 0
  outimg[!not_all_same, ] = sub_replacement
  stopifnot(!any(is.na(outimg)))


  ties.method = match.arg(ties.method)
  label = umat[max.col(outimg, ties.method = ties.method)]

  colnames(p) = colnames(q) = umat

  L = list(
    sensitivity = p,
    specificity = q,
    probability = outimg,
    label = label,
    prior = prior,
    number_iterations = iiter,
    convergence_threshold = tol,
    convergence_value = diff
  )
  return(L)
}
