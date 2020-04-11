## Identifies gene modules based on a simple mixture of von-Mises distributions
## using the intution that correlations are simply inner products of unit vectors

library('WGCNA')  # only used for labels2colors
library('movMF')
library('impute')
library('argparse')
library('parallel')

getArgs <- function() {
  parser <- ArgumentParser()
  parser$add_argument('expr', help='Expression matrix. Genes are rows, samples are columns.')
  parser$add_argument('out', help=paste('The output gene clustering file.',
                                        'Provides cluster assignment and cluster probability scores.'))
  parser$add_argument('--minK', help='Minimum number of clusters', default=5)
  parser$add_argument('--maxK', help='Maximum number of clusters', default=20)
  parser$add_argument('--transform', help='Data transformation to use: none, spearman, or bicor',
                      default='none')
  parser$add_argument('--cores', help='Number of cores', default=1, type='integer')
  parser$add_argument('--vmf_out', help='Output file for VMF models', default="NONE")
  parser$add_argument('--bic_mult', help='BIC multiplier', default=1.0, type='double')
  parser$add_argument('--no_bg', help='No background fit for gene classification', default=FALSE, type='logical')
  args <- parser$parse_args()
  
  if ( ! args$transform %in% c('spearman', 'bicor', 'none') ) {
    stop('Transform must be one of "spearman", "bicor", or "none"')
  }

  args
}

logUniformSphere <- function(n.dim) {
  (n.dim/2) * log(2*pi) - lgamma(n.dim/2)
}

scaledUnit <- function(v, center=T) {
  if ( center ) {
    v <- v - mean(v, na.rm=T)
  }
  v / sqrt(sum(v^2, na.rm=T))
}

spearman.transform <- rank

bicor.transform <- function(vec) {
  # transform the vector so that inner products will produce
  # the bi-weight midcorrelation
  m <- median(vec, na.rm=T)
  v <- mad(vec, na.rm=T)
  vec_qn <- (vec - m)/(9 * v)
  wgt <- (1 - vec_qn^2)^2 * ((1 - abs(vec_qn)) > 0)
  (vec - m) * wgt
}

loadExpression <- function(expr.file, transform) {
  print(expr.file)
  expr <- read.table(expr.file, header=T)
  e.sd <- apply(expr, 1, mad, na.rm=T)
  if ( any(e.sd <= 1e-6) ) {
    print(paste('Removing', sum(e.sd <= 1e-6), 'low-variance genes'))
    expr <- expr[e.sd > 1e-6,]
  }
  if ( transform == 'spearman' ) {
    expr <- t(apply(expr, 1, spearman.transform))
  } else if ( transform == "bicor" ) {
    expr <- t(apply(expr, 1, bicor.transform))
  }
  print(dim(expr))
  if ( rownames(expr)[1] == "1" ) {
    # no row-names. Try to assign them.
    if ( is.numeric(expr[1,1]) ) {
      stop('Expression data should have row-names, or the first column should consist of gene names.')
    }
    
    rownames(expr) <- expr[,1]
    expr <- expr[,-1]  # drop the names column
  }
  # now normalize the vectors to mean-0 unit-norm
  norm.expr <- t(apply(expr, 1, scaledUnit, center=(transform != 'bicor')))
  if ( any(is.na(norm.expr)) ) {
    print('NA values in expression. Imputing and re-normalizing.')
    messages <- capture.output({
      expr <- impute.knn(norm.expr)$data
    }, file=NULL)
    norm.expr <- t(apply(expr, 1, scaledUnit, center=(transform != 'bicor')))
  }
  norm.expr
}

addBackground <- function(data, vmf.fit) {
  # add a background uniform component to the data
  vmf.theta <- vmf.fit$theta  # mixture vMF parameters
  vmf.alpha <- vmf.fit$alpha  # mixture proportions
  new.theta <- rbind(vmf.theta, rep(1e-6, ncol(vmf.theta)))
  bkg.lik <- function(mix.mix.alpha) {
    new.alpha <- c((1 - mix.mix.alpha) * vmf.alpha, mix.mix.alpha)
    sum(dmovMF(data, new.theta, new.alpha, log=T)) + dbeta(mix.mix.alpha, 1, 2.5, log=T)
  }
  opt.res <- optimize(bkg.lik, lower=0.005, upper=0.995, maximum=T)
  print(paste('Unform background component fit at', opt.res[['maximum']]))
  new.alpha <- c(vmf.alpha * (1 - opt.res[['maximum']]), opt.res[['maximum']])
  new.lik <- opt.res[['objective']]
  # match what vMF does
  attr(new.lik, 'df') <- prod(dim(data)) - 1L - length(new.alpha)
  attr(new.lik, 'nobs') <- nrow(data)
  class(new.lik) <- 'logLik'
  res <- list(theta=new.theta, alpha=new.alpha, lik=new.lik,
              L=new.lik, P=vmf.fit$P, iter=vmf.fit$iter)
  class(res) <- 'movMF'
  res
} 

fitPredictVMF <- function(data, min.k, max.k, verbose=T, cores=1, add_bg=TRUE, bic.mul=1.0) {
  print(sprintf('Clustering %d points', nrow(data)))
  vMF.runs <- mclapply(min.k:max.k, function(k) {
    mixture <- movMF(data, k, reltol=1e-6, maxiter=1000, nruns=2, verbose=verbose) #minalpha=0.01) -- buggy // if last run dropped some components
    if ( add_bg ) {
      mixture.g <- addBackground(data, mixture)
      cluster.loglik <- mixture.g[['lik']]
    } else {
      mixture.g <- mixture
      cluster.loglik <- sum(dmovMF(data, mixture$theta, mixture$alpha, log=T))
    }
    mix.bic <- cluster.loglik - 2 * bic.mul * ncol(data) * (k + 1) * log(nrow(data))
    list(bic=mix.bic, model=mixture.g)
  }, mc.cores=cores)
  BIC <- lapply(vMF.runs, function(x) { x[['bic']] })
  mix.models <- lapply(vMF.runs, function(x) { x[['model']] })
  k.df <- data.frame(BIC=unlist(BIC), 
                     actual.k=sapply(mix.models, function(x) { length(x$alpha) - 1}),
                     initial.k=min.k:max.k)
  if ( verbose ) {
    print('Runs:')
    print(k.df)
    print("Best")
    print(k.df[which.max(k.df$BIC),])
  }
  best.bic <- which.max(unlist(BIC))
  best.mix <- mix.models[[best.bic]]
  best.mix.bic <- mix.models[[best.bic]][['bic']]
  best.k <- length(best.mix$alpha) - 1
  # make predictions
  hard.preds <- predict(best.mix, type='class_ids')
  # note: the unIform background is the last component
  #hard.preds[hard.preds == nrow(best.mix$theta)] <- 0
  soft.preds <- predict(best.mix, type='memberships')
  # remove the "Scattered" background
  #scatter.thresh <- exp(logUniformSphere(ncol(data)))
  #scattered.objects <- which(apply(soft.preds, 1, max) < scatter.thresh)
  #hard.preds[scattered.objects] <- 0
  colnames(soft.preds) <- sapply(colnames(soft.preds), function(z) {
    paste('mem', z, sep='.')
  })
  list('classes'     = hard.preds, 
       'memberships' = soft.preds, 
       'BIC'         = best.mix.bic,
       'mixture'     = best.mix,
       'k'           = best.k)
}

main <- function(args) {
  print(args)
  if ( ! args$no_bg ) {
    stop('Supply --no_bg TRUE')
  }
  data <- loadExpression(args$expr, args$transform)
  vMF.results <- fitPredictVMF(data, args$minK, args$maxK, verbose=T, args$cores, ! args$no_bg, args$bic_mul)
  out.data <- data.frame(vMF.results[['memberships']])
  out.data$module <- vMF.results[['classes']]
  rownames(out.data) <- rownames(data)
  write.table(format(out.data, scientific=T, digits=4), file=args$out, quote=F)
  if ( args$vmf_out != 'NONE' ) {
    model.params <- vMF.results[['mixture']]$theta
    print(names(vMF.results[['mixture']]))
    print(dim(model.params))
    colnames(model.params) <- colnames(data)
    m.alpha <- as.matrix(vMF.results[['mixture']]$alpha)
    colnames(m.alpha) <- 'alpha'
    model.params <- cbind(m.alpha, model.params)
    rownames(model.params) <- labels2colors(c(1:vMF.results[['k']], 0))
    write.table(format(model.params, scientific=F, digits=4), file=args$vmf_out, quote=F)
  } 
}


if ( ! interactive() ) {
  args <- getArgs()
  main(args)
}
