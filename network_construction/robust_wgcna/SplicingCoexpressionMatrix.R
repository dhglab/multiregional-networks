# compute co-expression based on splicing, using
# principal angles. 
require('expm')
require('parallel')
library('data.table')

vec.bicor.transform <- function(vec) {
  # transform the vector so that inner products will produce
  # the bi-weight midcorrelation
  m <- median(vec, na.rm=T)
  v <- mad(vec, na.rm=T)
  vec_qn <- (vec - m)/(9 * v)
  wgt <- (1 - vec_qn^2)^2 * ((1 - abs(vec_qn)) > 0)
  q <- (vec - m) * wgt
}

vec.spearman.transform <- frank

mat.bicor.transform <- function(mat, center=T) {
  M = matrix(vec.bicor.transform(as.vector(mat)), nrow=nrow(mat), byrow=F)
  if ( center ) {
    M = t(scale(t(M), scale=F))
  }
  as.matrix(M)
}

mat.spearman.transform <- function(mat, center=T) {
  M = matrix(vec.spearman.transform(as.vector(mat)), nrow=nrow(mat), byrow=F)
  if ( center ) {
    M = t(scale(t(M), scale=F))
  }
  as.matrix(M)
}

bicor.transform <- function(expr, tx.map, center=T, threads=1) {
  tx.map.sub <- subset(tx.map, isoform %in% rownames(expr))
  genes <- unique(tx.map.sub$gene)
  if ( threads == 1 ) {
    for ( gene in genes ) {
      trans <- subset(tx.map.sub, gene == gene)$isoform
      tx <- mat.bicor.transform(expr[trans,], center=center)
      expr[trans,] <- tx
    }
  } else {
    iso.list <- mclapply(genes, function(gene) {
      trans <- subset(tx.map.sub, gene == gene)$isoform
      mat.bicor.transform(expr[trans,], center=center)
    })
    expr <- do.call(rbind, iso.list)
  }
  expr
}

spearman.transform <- function(expr, tx.map, center=T, threads=1) {
  tx.map.sub <- subset(tx.map, isoform %in% rownames(expr))
  genes <- unique(tx.map.sub$gene)
  if ( threads == 1 ) {
    for ( gene in genes ) {
      trans <- subset(tx.map.sub, gene == gene)$isoform
      tx <- mat.spearman.transform(expr[trans,], center=center)
      expr[trans,] <- tx
    }
  } else {
    iso.list <- mclapply(genes, function(gene) {
      trans <- subset(tx.map.sub, gene == gene)$isoform
      mat.spearman.transform(expr[trans,], center=center)
    })
    expr <- do.call(rbind, iso.list)
  }
  expr
}

cond.number <- function(M, tol=1e-6, steps=50) {
  # approximate the condition number of M
  q <- rnorm(nrow(M))
  q <- q/sqrt(sum(q^2))
  v1 <- v2 <- q
  for ( s in 1:steps ) {
    v1.t <- M %*% v1
    r1 <- t(v1) %*% v1.t
    v2.t <- M %*% v2
    r2 <- t(v2) %*% v2.t
    v1.n <- v1.t - r1 * v1
    v2.n <- r2 * v2 - v2.t
    v1.n <- v1.n/sqrt(sum(v1.n^2))
    v2.n <- v2.n/sqrt(sum(v2.n^2))
    if ( sqrt(sum((v1.n-v1)^2)) + sqrt(sum((v2.n-v2)^2)) < tol ) {
      break
    }
    v1 <- v1.n
    v2 <- v2.n
  }
  r1/r2
}

inv.sqrtm <- function(A) {
  # computes the inverse square-root of A
  # (faster than SVD method for large A)
  if ( nrow(A) == 1 ) {
    1/sqrt(A)
  } else {
    A.log <- logm(A, method='Higham08')
    expm(-0.5 * A.log)
  }
}

semi.orth <- function(mat) {
  # convert `mat` (k x n with n > k) to a 
  # semi-orthogonal matrix
  # [ mat %*% t(mat) == I ]
  MM <- mat %*% t(mat)
  # probably should check for singularity
  if ( nrow(MM) > 1 && kappa(MM, exact=F) > 1e6 ) {
    d.o <- diag(MM)
    MM.c <- cov2cor(MM) + diag(rep(1e-9, nrow(MM)))
    MM.c <- cov2cor(MM.c)
    for ( i in 1:nrow(MM.c) ) {
      for ( j in 1:ncol(MM.c) ) {
        MM[i, j] <- sqrt(d.o[i] * d.o[j]) * MM.c[i,j]
      }
    }
  }
  inv.sqrtm(MM) %*% mat
}

softAdj <- function(cos.angles, power=6) {
  # extension of soft adjacency to hyperplanes
  # this looks arbitrarily like the geometric mean
  # of each soft adjacency, but it is the same
  # soft transformation applied to the determiant
  # of the "hyperplane angle" (i.e. extension of
  # inner product)
  # 
  # to be fair there are two definitions:
  #  (0.5 + 0.5 * det(A))^p
  # and
  # det([diag(0.5) + diag(0.5) * A])^p
  # the former 
  sft = (0.5 + 0.5 * cos.angles)^power
  prod(sft)^(1/length(sft)) 
}

orthoNormalize <- function(iso.expr, gene.iso.map) {
  # make every gene's isoformd data semi-orthogonal
  iso.orth <- iso.expr
  for ( gene in unique(gene.iso.map$gene) ) {
    isos <- gene.iso.map$isoform[gene.iso.map$gene == gene]
    iso.orth[isos,] <- semi.orth(iso.expr[isos,,drop=F])
    check <- iso.orth[isos,,drop=F] %*% t(iso.orth[isos,,drop=F])
    if ( ! all(abs(diag(check) - 1) < 1e-4) ) {
      print(gene)
      print(check)
      print(diag(check) - 1)
      stop('Did not normalize...')
    }
  }
  iso.orth
}

checkIBAInputs <- function(isoform.expr, gene.isoform.map, method) {
  if ( ! method %in% c('pearson', 'spearman', 'bicor', 'raw') ) {
    stop('Method must be one of "raw", "pearson", "spearman", or "bicor"')
  }

  if ( method == 'spearman' ) {
    print('Warning: Spearman transform will take some time')
    isoform.expr <- rank.transform(isoform.expr, gene.isoform.map)
  } else if ( method == 'bicor' ) {
    isoform.expr <- bicor.transform(isoform.expr, gene.isoform.map)
  } else if ( method == 'pearson' ) {
    isoform.expr <- t(scale(t(isoform.expr)))
  } else {
    isoform.expr <- t(scale(t(isoform.expr), scale=F))
  }

  print('Checking isoform data for singular quantifications...')

  drop.isos <- list()
  for ( gene in unique(gene.isoform.map$gene) ) {
    isos <- gene.isoform.map$isoform[gene.isoform.map$gene == gene]
    check <- isoform.expr[isos,,drop=F] %*% t(isoform.expr[isos,,drop=F])
    this.drop <- list()
    iso.means <- apply(isoform.expr[isos,,drop=F], 1, mean)
    while ( kappa(check, exact=F) > 1e6 ) {
      # singular matrix, start dropping isoforms by expression
      least.expr <- which.min(iso.means)
      this.drop[[1 + length(this.drop)]] <- rownames(check)[least.expr]
      isos <- setdiff(isos, unlist(this.drop))
      check <- isoform.expr[isos,,drop=F] %*% t(isoform.expr[isos,,drop=F])
      iso.means <- apply(isoform.expr[isos,,drop=F], 1, mean)
    }
    if ( length(this.drop) > 0) {
      drop.isos <- c(drop.isos, this.drop)
      print(sprintf('Dropped %d for %s final size %d', length(this.drop), gene, nrow(check)))
    }
  }
  which.drop <- which(rownames(isoform.expr) %in% unlist(drop.isos))
  isoform.expr[-which.drop,,drop=F]
}


isoformBasedAdjacencyParallel <- function(isoform.expr, gene.isoform.map, method='raw', verbose=T, ncores=2, njobs=40) {
  if ( ncores < 2 ) {
    stop('Must have more than 1 core')
  }
  isoform.expr <- checkIBAInputs(isoform.expr, gene.isoform.map, method)
  isoform.expr <- orthoNormalize(isoform.expr, gene.isoform.map)
  genes <- unique(gene.isoform.map$gene)
  if ( verbose ) {
    print('Pre-computing indeces')
  }
  iso.idx <- lapply(genes, function(gene) { 
    which(rownames(isoform.expr) %in% gene.isoform.map$isoform[gene.isoform.map$gene == gene])
  })
  names(iso.idx) <- genes
  # cut the matrix calculation into chunks. There are
  #   n * (n-1)/2 pairs of the form
  #   1 2
  #   1 3
  #   ...
  #   1 n
  #   2 3
  #   2 4
  #   ...
  #   2 n
  #   ...
  # n-1 n
  #
  #  first define a function that maps the kth row
  #  to the gene idx it represents
  idxToCoords <- function(k, n) {
    # let the first coordinate be `r`
    # there are `(n-1) + (n-2) + ... + (n - r)` rows up to
    # in which the first coordinate is `r` or smaller
    # this sum is
    #   0.5 * r * (2 * n - r - 1) = C_r
    # choosing the largest `r` such that C_r < k gives
    # the first coordinate. The solution at equality is
    #   r = 0.5 * ( 2 * n - sqrt((1 - 2*n)^2 - 8*k) - 1 )
    # and so
    rad = sqrt((1 - 2 * n)^2 - 8 * k)
    r = as.integer(ceiling(0.5 * (2 * n - rad - 1)))
    # then the number of additional rows in determines
    # the column; as we will have
    #  c = 1 + r + `(k - C_r)`
    C_r = 0.5 * (r-1) * (2 * n - (r-1) - 1)
    cl = r + (k - C_r)
    c(r, cl)
  }
  # next calculate the length of the vectorized triangular matrix
  tlen = length(genes) * (length(genes) - 1)/2
  # cut it into `njobs` jobs
  job.idx <- as.integer(seq(1, tlen, length.out=njobs))
  # map each job to a submatrix
  if ( verbose ) {
    print(sprintf('Computing pairwise Stiefel distance using %d cores', ncores))
  }
  iba.matrix <- mclapply(2:njobs, function(jnum) {
    print(sprintf('job %d/%d..', jnum, njobs))
    jidx.start <- job.idx[jnum-1]
    jidx.end <- job.idx[jnum] - (jnum < njobs)  # clopen except last interval
    jcoords.start <- idxToCoords(jidx.start, length(genes))
    jcoords.end <- idxToCoords(jidx.end, length(genes))
    rc.vec <- rep(0, jidx.end-jidx.start+1)
    # we know that these coordinates are sequential
    row <- jcoords.start[1]
    col <- jcoords.start[2]
    row.mat <- col.mat <- NULL
    jidx <- 1
    print(sprintf('@@ (j=%d, outer loop)', jnum))
    while ( row <= jcoords.end[1] ) {
      row.gene <- genes[row]
      row.iso <- iso.idx[[row.gene]]
      print(sprintf('@@ (j=%d, inner loop)', jnum))
      row.mat <- isoform.expr[row.iso,,drop=F]
      while ( col <= length(genes) ) {
        col.gene <- genes[col]
        col.iso <- iso.idx[[col.gene]]
        col.mat <- isoform.expr[col.iso,,drop=F]
        rc.vec[jidx] <- softAdj(svd(row.mat %*% t(col.mat))$d)
        col <- col + 1
        jidx <- jidx + 1
        if ( jidx > length(rc.vec) ) {
          break
        }
        print(sprintf('@@ (j=%d, inner loop col %d < %d, jidx %d <= %d)', jnum, col, length(genes), jidx, length(rc.vec)))
      }
      print(sprintf('@@ (j=%d, inner loop is done)', jnum))
      if ( jidx > length(rc.vec) ) {
        break
      }
      row <- row + 1
      col <- row + 1
      print(sprintf('@@ (j=%d r=%d/%d c=%d)', jnum, row, jcoords.end[1], col))
    }
    print(sprintf('job %d done...', jnum))
    rc.vec
  }, mc.cores=ncores, mc.preschedule = FALSE)
  A <- matrix(0, nrow=length(genes), ncol=length(genes))
  rownames(A) <- colnames(A) <- genes
  A[lower.tri(A, diag=F)] <- unlist(iba.matrix)
  A <- t(A)
  A <- A + t(A)
  diag(A) <- 1
  A
}

isoformBasedAdjacency <- function(isoform.expr, gene.isoform.map, method='raw', verbose=T, ncores=1, pow=6, njobs=40) {
  if ( ncores > 1 ) {
    return(isoformBasedAdjacencyParallel(isoform.expr, gene.isoform.map, method, verbose, ncores, njobs=njobs))
  }
  isoform.expr <- checkIBAInputs(isoform.expr, gene.isoform.map, method)
  gene.isoform.map <- subset(gene.isoform.map, isoform %in% rownames(isoform.expr))
  if ( verbose ) {
    print('Normalizing to identity')
  }
  isoform.expr <- orthoNormalize(isoform.expr, gene.isoform.map)
  # pre-compute isoform indeces
  genes <- unique(gene.isoform.map$gene)
  if ( verbose ) {
    print('Pre-computing indeces')
  }
  iso.idx <- lapply(genes, function(gene) { 
    which(rownames(isoform.expr) %in% gene.isoform.map$isoform[gene.isoform.map$gene == gene])
  })
  names(iso.idx) <- genes
  gene.adj <- matrix(0, nrow=length(genes), ncol=length(genes))
  rownames(gene.adj) <- colnames(gene.adj) <- genes
  if ( verbose ) {
    print('Computing pairwise Stiefel distance')
  }
  for ( ix1 in 1:(length(genes)-1) ) {
    G1 <- isoform.expr[iso.idx[[ix1]],,drop=F]
    for ( ix2 in (1+ix1):length(genes) ) {
      G2 <- isoform.expr[iso.idx[[ix2]],,drop=F]
      # ith element = cosine of ith principal angle
      principal.cors <- svd(G1 %*% t(G2), nu=0, nv=0)$d
      if ( any(abs(principal.cors) > 1.1) ) {
        stop(sprintf('bad correlation for %s-%s', genes[ix1], genes[ix2]))
      } 
      if ( any(abs(principal.cors) > 1) ) {
        principal.cors[principal.cors > 1] <- 1
        principal.cors[principal.cors < -1] <- -1
      }
      #gene.dist[ix1, ix2] <- gene.dist[ix2,ix1] <- sqrt(sum(acos(principal.angles)^2))
      gene.adj[ix1, ix2] <- gene.adj[ix2, ix1] <- softAdj(principal.cors, power=pow)
    }
    if ( (ix1 %% 50 == 0) &&  verbose ) {
      k <- ix1 * (ix1-1)/2
      n <- length(genes) * (length(genes)-1)/2
      print(sprintf('Done: %d/%d genes (%d/%d pairs)', ix1, length(genes), k, n))
    }
  }
  diag(gene.adj) <- 1
  gene.adj
}

