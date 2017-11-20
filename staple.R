##############################
#
##############################
rm(list=ls())
library(matrixStats)
library(RNifti)

ensure_Nifti = function(x) {
	if (is.list(x)) {
		x = lapply(x, ensure_Nifti)
	}
	if (is.character(x)) {
		if (length(x) > 0) {
			x = lapply(x, readNifti)
		} else {
			x = readNifti(x)
		}
	}
	x
}

files = list.files(
	pattern = "Manual.*.nii.gz")
# imgs = check_nifti(files)
imgs = ensure_Nifti(files)
orientations = sapply(imgs, orientation)
ori = orientation(imgs[[1]])

imgs = lapply(imgs,
	function(x) {
		orientation(x) = ori
		x
	})

orig = t(sapply(imgs, c))

# bad = which(orig > 1, arr.ind = TRUE)

orig = orig > 0
stopifnot(!any(is.na(orig)))
keep = colSums(orig)
stopifnot(!any(is.na(keep)))


mat = orig[, keep > 0]
umat = unique(c(mat))
f_t_i = colMeans(mat, na.rm = TRUE)

dmat = 1 - mat

mat[mat==0] = NA
dmat[dmat==0] = NA


n_readers = nrow(mat)
n_voxels = ncol(mat)
###################
#initialize
p = q = rep(0.99999, n_readers)

max_iter = 1000
tol = .Machine$double.eps


# mat is D
### run E Step
for (i in seq(max_iter)) {
	pmat = p * mat 
	pmat = colProds(pmat, na.rm = TRUE)
	sep_pmat = (1-p) * dmat
	sep_pmat = colProds(sep_pmat, na.rm = TRUE)

	qmat = q * mat
	qmat = colProds(qmat, na.rm = TRUE)
	sep_qmat = (1-q) * dmat
	sep_qmat = colProds(sep_qmat, na.rm = TRUE)

	a_i = f_t_i * pmat * sep_pmat
	b_i = (1-f_t_i) * qmat * sep_qmat
	W_i = a_i/(a_i + b_i)


	sum_w = sum(W_i)

	new_p  = t(mat) * W_i
	new_p  = colSums(new_p,	na.rm = TRUE)
	new_p = new_p/sum_w

	new_q  = t(dmat) * (1 - W_i)
	new_q  = colSums(new_q,	na.rm = TRUE)
	new_q = new_q/(n_voxels - sum_w)

	diff_p = abs(p - new_p)
	diff_q = abs(q - new_q)
	diff = max(c(diff_p, diff_q))
	if (diff <= tol) {
		print("Convergence!")
		break
	} else {
		print(paste0("iter: ", i, 
			", diff: ", diff))
	}

	p = new_p
	q = new_q 
}

stopifnot(!any(is.na(W_i)))

outimg = rep(0, ncol(orig))
outimg[keep > 0] = W_i
outimg = array(outimg, 
	dim = dim(imgs[[1]]))

hdr = RNifti::dumpNifti(imgs[[1]])
hdr$cal_max = 1
hdr$cal_min = 0
hdr$dattype = 16

outimg = RNifti::updateNifti(
	outimg, template = hdr)



thresh = outimg >= 0.5
center = neurobase::xyz(thresh)

flair = readNifti("3DFLAIR.nii.gz")
orientation(flair) = ori

neurobase::ortho2(flair, outimg, xyz = center)

