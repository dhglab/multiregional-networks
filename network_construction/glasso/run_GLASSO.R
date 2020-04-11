library('Matrix')
library('parallel')
library('QUIC')
library('WGCNA')

args <- commandArgs(trailingOnly=T)
if ( length(args) < 3 ) {
  stop('Usage: run_GLASSO.R [data.file] [output.Rda file] [num.cores]')
}
data.file <- args[1]
output.file <- args[2]
n.cores <- as.integer(args[3])

print(sprintf('Will save to %s', output.file))

run.blockwise.GLASSO <- function(expr.data, block.size=400, n.boot=20, beta=0.05, seed=20160428, n.cores=1, max.iter=20) {
  ## :expr.data:    data matrix; rows are samples, columns are genes, and we want to cluster columns
  ## :block.size:   use PAM to pre-cluster covariance matrix into blocks of (approximately) this size
  ## :n.boot:       number of bootstraps for instability estimation
  ## :beta:         target instability for lambda parameter
  ## :seed:         seed for bootstrapping
  ## :n.cores:      number of cores to use
  ## :max.iter:     maximum iterations for refining the lambda parameter
  #####
  set.seed(seed)
  if ( n.cores > 1 ) {
    enableWGCNAThreads(n.cores)
  }
  # expr.data :: columns are genes and rows are samples
  # pre-cluster into blocks using projective k-means (WGCNA)
  inc = as.integer(ncol(expr.data)/block.size)*3
  clustering = projectiveKMeans(expr.data, preferredSize=block.size, checkData=T, sizePenaltyPower=5,
                                nCenters=inc, verbose=5)
  if ( n.cores > 1 ) {
    disableWGCNAThreads()
  }

  # mcable
  res <- lapply(unique(clustering$clusters), function(block) {
    init.lambda <- find.init.lambda(expr.data[,clustering$clusters==block], init.lambda=25, name=sprintf('block %s', block))
    sub.size <- as.integer(0.8 * nrow(expr.data))
    stars.lambda <- run.GLASSO.STARS(expr.data[,clustering$clusters==1], init.lambda, sub.size, n.boot, beta)
    gcov <- cov(expr.data[,clustering$clusters==block])
    print(paste('clustering block', block, 'cov=', dim(gcov), 'lambda=', stars.lambda['Lambda']))
    #glasso(gcov, 1/stars.lambda['Lambda'])
    list('QUIC'=QUIC(gcov, rho=1/stars.lambda[['Lambda']]), 'index'=which(clustering$clusters == block),
         'block'=block)
  })#, mc.cores=n.cores)
  res
}



find.init.lambda <- function(expr.data, min.prop=0.05, max.prop=0.2, init.lambda=25, name='') {
  # finds an initial value of lambda that 
  #print('finding initial value for Lambda')
  gcov <- cov(expr.data)
  Lambda <- init.lambda
  dL <- 0.8 * Lambda
  #gl.res <- glasso(gcov, 1/Lambda, approx=T)
  gl.res <- QUIC(gcov, rho=1/Lambda)
  nc <- ncol(expr.data)
  p.edges <- sum(abs(as.dist(gl.res[['W']]) > 1e-10))/(nc * (nc - 1)/2)
  #print(paste('starting: lambda ', Lambda, 'p.edges', p.edges))
  while ( p.edges > max.prop ) {
    print(paste(name, ' - Step1: Lambda=',Lambda,'p.edges=',p.edges))
    gl.res <- QUIC(gcov, rho=1/Lambda)
    p.edges <- sum(abs(as.dist(gl.res[['W']])) > 1e-10)/(nc * (nc - 1)/2)
    Lambda = Lambda - dL
    dL <- 0.8 * Lambda
  }
  
  while ( p.edges < min.prop ) {
    print(paste(name, ' - Step2: Lambda=',Lambda,'p.edges=',p.edges))
    Lambda = Lambda + dL
    gl.res <- QUIC(gcov, rho=1/Lambda)
    #p.edges <- sum(1 * (abs(gl.res[['W']]) > 1e-10))/(ncol(expr.data)^2)
    p.edges <- sum(abs(as.dist(gl.res[['W']])) > 1e-10)/(nc * (nc - 1)/2)
  }

  dL = Lambda/10
  #print(sprintf('dL=%f', dL))
  
  iters = 1 
  previous = '?'
  prev.Lambda <- -1
  while ( (p.edges < min.prop || p.edges > max.prop) && iters < 20 ) {
    print(paste(name, ' - Step3: Lambda=',Lambda,'p.edges=',p.edges))
    if ( p.edges < min.prop ) {
      if ( previous == '+' ) {
        # take average
        h <- Lambda
        Lambda <- (prev.Lambda + Lambda)/2
        prev.Lambda <- h
        dL <- abs(Lambda - prev.Lambda)
      } else {
        prev.Lambda <- Lambda
        Lambda <- Lambda + dL
      }
      previous = '-'
    } else {
      if ( previous == '-' ) {
        h <- Lambda
        Lambda <- (prev.Lambda + Lambda)/2
        prev.Lambda <- h
        dL <- abs(Lambda - prev.Lambda)
      } else {
        prev.Lambda <- Lambda 
        Lambda <- Lambda - dL
      }
      previous = '+'
    }
    if ( Lambda <= 0 ) {
      save(gcov, file='init.debug.Rda')
      stop('Initial value of lambda could not be found. See init.debug.Rda')
    }
    gl.res <- QUIC(gcov, rho=1/Lambda)
    #p.edges <- sum(1 * (abs(gl.res[['W']]) > 1e-10))/(ncol(expr.data)^2)
    p.edges <- sum(abs(as.dist(gl.res[['W']])) > 1e-10)/(nc * (nc - 1)/2)
    iters = 1 + iters
  }

  print(sprintf('%s Init:: Lambda=%f, p.edges=%f', name, Lambda, p.edges))
  Lambda
}

run.STARS.step <- function(data, boot.seeds, sub.size, Lambda, n.cores=1) {
  sample.sets <- lapply(boot.seeds, function(seed) {
    sample.int(dim(data)[1], size=sub.size)
  })
  ec.agg <- NULL
  start <- 1
  end <- start + n.cores - 1
  while ( start <= length(sample.sets) ) {
    emax <- min(end, length(sample.sets))
    # mcable
    edge.counts <- mclapply(start:emax, function(idx) {
      run.GLASSO.boot(data, sample.sets[[idx]], sub.size, Lambda)
    }, mc.cores=n.cores, mc.preschedule=F)
    if ( is.null(ec.agg) ) {
      ec.agg <- Reduce('+', edge.counts)/length(boot.seeds)
      rm(edge.counts)
    } else {
      ec.agg <- ec.agg + Reduce('+', edge.counts)/length(boot.seeds)
      rm(edge.counts)
    }
    start <- start + n.cores
    end <- start + n.cores - 1
  }
  print(dim(ec.agg))
  vars <- 2 * ec.agg * (1 - ec.agg)
  print(sprintf('Lambda=%f, mean edge weight=%f', Lambda, mean(as.dist(ec.agg))))
  instab <- sum(vars[upper.tri(vars)])/(ncol(data)*(ncol(data)-1)/2)
  instab
}


run.GLASSO.STARS <- function(expr.data, init.lambda=NULL, sub.size=40, n.boot=100, beta=0.1) {
  # run GLASSO algorithm with STARS selection for the penalty term, see
  # Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models
  # (Liu H, Roeder K, Wasserman L; NIPS 2010)
  if ( is.null(init.lambda) ) {
    dfactor <- sum(cov(expr.data))
    Lambda <- (1/beta) * sqrt(dfactor)
  } else {
    Lambda <- init.lambda
  }
  print(paste('Running STARS with initial value of ', init.lambda))
  dL <- Lambda/50
  Lambda.grid <- seq(0.5*Lambda, 1.2*Lambda, by=dL)
  init.W <- QUIC(cov(expr.data), rho=1/Lambda)[['W']]
  edges <- lapply(Lambda.grid, function(lg) {
    rep(0, length(as.dist(init.W)))
  })
  for ( boot in 1:n.boot ) {
    boot.cov <- cov(expr.data[sample.int(nrow(expr.data), sub.size),])
    edge.grid <- lapply(Lambda.grid, function(lmb) {
      gres <- QUIC(boot.cov, rho=1/lmb, W.init=init.W, msg=0)
      as.dist(1 * abs(gres[['W']]) > 1e-10)
    })#, mcores=n.cores)
    for ( i in 1:length(edges) ) {
      edges[[i]] = edges[[i]] + edge.grid[[i]]
    }
  }
  nc <- ncol(expr.data)
  instability <- sapply(edges, function(ecounts) {
      efrac = ecounts/n.boot
      evar = 2 * efrac * (1 - efrac)
      sum(evar)/(nc * (nc - 1) / 2)
  })

  mean.edge <- sapply(edges, function(ecounts) {
      efrac = ecounts/n.boot
      mean(efrac)
  })
  
  best.idx <- which.min(abs(instability - beta))
  Lambda.dat <- data.frame(Lambda=Lambda.grid, instability=instability, mean.edge=mean.edge)
  print(Lambda.dat)
  tr = list('Lambda'=Lambda.grid[best.idx], 'instab'=instability[best.idx])
  #print(tr)
  tr
}
                                                                      

run.GLASSO.boot <- function(data, samples, size, Lambda) {
  gcov <- cov(data[samples,])
  gres <- QUIC(gcov, rho=1/Lambda)
  1 * (abs(gres[['W']]) > 0)
}

spectral.clust <- function(Q, normalize=T) {
  diag(Q) <- 0
  Q[Q < 0] <- 0
  DQ <- rowSums(Q)
  background <- which(DQ == 0)

  if ( length(background) > 0 ) {
    QP <- Q[-background, -background]
    DQP <- DQ[-background]
  } else {
    QP <- Q
    DQP <- DQ
  }

  if ( normalize ) {
    k <- length(DQP)
    # symmetric normalized laplacian
    L <- diag(rep(1, k)) - diag(1/sqrt(DQP)) %*% QP %*% diag(1/sqrt(DQP))
  } else {
    L <- DQP - QP
  }

  if ( all(L < 1e-14 & L > -1e-14) ) {
    return(rep(0, nrow(Q)))
  }

  L.svd <- eigen(L, symmetric=T)
  null.cutoff <- 1e-12
  n.clusters <- sum(L.svd$values < null.cutoff)
  if ( n.clusters > 1 ) {
    rep.vec.idx <- (nrow(L):1)[(1+n.clusters):(2*n.clusters)]
    rep.vec <- L.svd$vectors[,rep.vec.idx]
  
    clusters <- kmeans(rep.vec, n.clusters)$cluster
  } else {
    clusters <- rep(1, nrow(L))
  }
  clusters.final <- rep(0, nrow(Q))
  if ( length(background) > 0 ) {
    clusters.final[-background] <- clusters
    clusters.final[background] <- 0
  } else {
    clusters.final <- clusters
  }
  clusters.final
}

soft.cluster <- function(data, hard.clustering, dist.quantile=0.1, conf.t=0.75, A=NULL) {
  ## :input data:             raw expression data
  ## :input hard.clustering:  assignment of columns to clusters, 
  ##                          with 0 as unassigned
  ## :input dist.quantile:    The quantile of distance to use for epsilon-nn, will retain
  ##                          this proportion of edges
  ## :input conf.t:           The confidence threshold (proportion of NN for assignment)
  ## :input A:                Precomputed NN-graph
  if ( is.null(A) ) {
    d.mat <- dist(data, 'euclidean')
    d.thresh <- quantile(as.vector(d.mat), dist.quantile)
    A <- as.matrix(d.mat) <= d.thresh
    diag(A) <- F
  }
  new.clustering <- hard.clustering
  for ( z.idx in which(hard.clustering == 0) ) {
    nn.clust <- hard.clustering[which(A[z.idx,])]
    nn.clust <- nn.clust[nn.clust != 0]
    if ( length(nn.clust) == 0 ) {
      next
    }
    clust.count <- table(nn.clust)
    best.clust <- which.max(clust.count)
    confidence <- clust.count[best.clust]/sum(clust.count)
    if ( confidence > conf.t ) {
      new.clustering[z.idx] <- names(clust.count)[best.clust]
    }
  }
  new.clustering

}


cluster.GLASSO <- function(data, block.size=1200, n.boot=25, beta=0.05, 
                           seed=20160428, n.cores=1, max.iter=20,
                           dist.quantile=0.1, conf.t=0.75) {
  glasso.blocks <- run.blockwise.GLASSO(data, block.size, n.boot, beta, seed, n.cores, max.iter)
  print('Assigning clusters using spectral clustering + NN-propagation')
  block.clusters <- lapply(glasso.blocks, function(blk) {
    blk.pcov <- blk[[1]]$W
    blk.idx <- blk[[2]]
    blk.num <- blk[[3]]
    block.dat <- data[,blk.idx]  # columns genes, rows samples
    block.hard <- spectral.clust(blk.pcov)
    block.soft <- block.hard
    p.bs <- rep(-1, length(block.soft))
    while ( any(block.soft != p.bs) ) {
      p.bs <- block.soft
      block.soft <- soft.cluster(t(block.dat), block.soft, 
                                 dist.quantile, conf.t)
    }
    list('clusters'=sprintf('%d.%d', blk.num, block.hard), 
         'propa.clusters'=sprintf('%d.%s', blk.num, block.soft))
  })#, mcores=n.cores)
  print('Unwrapping blocks')
  block.idx <- do.call(c, lapply(glasso.blocks, function(b) b[[2]]))
  genes <- colnames(data)[block.idx]
  precision <- bdiag(lapply(glasso.blocks, function(b) b[[1]]$W))
  covariance <- bdiag(lapply(glasso.blocks, function(b) b[[1]]$X))
  hard.clust <- do.call(c, lapply(block.clusters, function(b) b[[1]]))
  soft.clust <- do.call(c, lapply(block.clusters, function(b) b[[2]]))
  save(list=ls(), file='debug.Rda')
  names(hard.clust) <- genes
  names(soft.clust) <- genes
  rownames(precision) <- colnames(precision) <- genes
  rownames(covariance) <- colnames(covariance) <- genes
  list('data.ordered'=data[,block.idx],
       'covariance'=covariance,
       'precision'=precision,
       'clusters'=hard.clust,
       'clusters.propagated'=soft.clust)
} 
  

data <- read.table(data.file, header=T, comment.char='')
if ( dim(data)[1] > dim(data)[2] ) {
  data <- t(data)  # columns: genes, rows: samples
}
data <- scale(data)

time.start <- proc.time()
#res <- run.blockwise.GLASSO(data, n.cores=n.cores)
res <- cluster.GLASSO(data, n.cores=n.cores)
time.end <- proc.time()

print('Start cluster')
print(time.start)
print('End cluster')
print(time.end)

print(names(res))

save(res, file=output.file)

