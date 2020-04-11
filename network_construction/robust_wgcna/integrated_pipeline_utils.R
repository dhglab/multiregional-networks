## utilities for the integrated pipeline
require('iterators')
require('itertools')

fstream <- function(filename, what, header=T, nlines=1) {
  handle = file(filename)
  open(handle)
  hdr <- NULL
  if ( header ) {
    hdr <- readLines(handle, n=1)
  }
  nxt <- function() {
    x <- scan(handle, what=what, nlines=nlines, quiet=T)
    if ( length(x) == 0 ) {
      stop('StopIteration')
    }
    x
  }
  obj <- list(nextElem=nxt)
  class(obj) <- c('fstream', 'abstractiter', 'iter')
  ihasNext(obj)
  #obj
}

## robust wgcna functions

unwrap.symm <- function(M) {
  # unwrap the symmetric matrix M into a vector consisting of the diagonal and upper triangular matrix
  d <- diag(M)
  ut <- as.vector(M[upper.tri(M)])
  c(d, ut)
}

sym.size.from.vec.size <- function(v) {
  # we have: z = len(v); k=dim(m)[1]; v = k + k*(k-1)/2 ;; k^2/2 + k/2 - v = 0
  # k = -(1/2) + sqrt(1/4 + 4*v/2)
  as.integer(round(sqrt(1/4 + 4*v/2) - 1/2))
}

wrap.symm <- function(v) {
  # wraps the vector v, consisting of the diagonal and upper triangular part of a symmetric matrix
  # back into a symmetric matrix
  k <- sym.size.from.vec.size(length(v))
  M <- matrix(0, nrow=k, ncol=k)
  diag(M) <- v[1:k]
  M[upper.tri(M)] <- v[-(1:k)]
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}
