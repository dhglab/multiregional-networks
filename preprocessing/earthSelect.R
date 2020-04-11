library('earth')
library('parallel')

earthSelect <- function(gene.expr, covars, n.predictors=100, num.genes=50, min.span=0, end.span=0, 
                        force.linear=T, nprune=NULL, max.terms=NULL, no.cross=NULL) {
  # uses the packages `earth` to determine an appropriate linear model for expression, given
  # technical covariates.
  # Inputs
  #  gene.expr     - The gene expression matrix (genes x samples)
  #  covars        - The sample covariate matrix (samples x covariates)
  #  n.predictors  - The (maximum) number of predictors to use in the forward earth  model
  #  num.genes     - The number of genes to sample for the model (max 1,000)
  #  min.span      - if force.linear=F, the number of points between any two spline knots
  #  end.span      - if force.linear=F, the number of points between a terminal knot and the end
  #                  of that particular covariate vector
  #  nprune        - Maximum number of (expanded) terms to keep in the backwards pass
  #  max.terms     - Maximum number of (expanded) terms to keep in the `final model`, e.g.
  #                  if there are 13 batches, but only batch2 and batch4 are selected
  #                  nprune counts this as 2 variables; while max.terms counts all 13
  #                  (since the -whole- factor needs to be included in the linear model)
  #  no.cross      - vector of covariates to exclude from cross terms
  #
  # Returns:
  #  earth         - the fitted `earth` object
  #  formula       - formula giving the right-hand-size of the covariate-correction LM
  #  terms         - the terms included
  #  term.imp      - importance for each term: this is the 80th percentile of max(beta)
  #  model.matrix  - the model matrix resulting from model.matrix(formula, covars)
  if ( is.null(max.terms) ) {
    max.terms <- 1000
  }
  if ( num.genes < nrow(gene.expr) ) { 
    gene.idx <- sample.int(nrow(gene.expr), size=num.genes)
  } else {
    gene.idx <- 1:nrow(gene.expr)
  }
  lhs <- paste('cbind(', paste(rownames(gene.expr)[gene.idx], collapse=', '), ')', sep='')
  rhs1 <- paste(colnames(covars), collapse=' + ')
  num.covs <- sapply(1:ncol(covars), function(c.idx) {
    ! is.factor(covars[,c.idx]) && ! is.character(covars[,c.idx])
  })
  num.covs <- colnames(covars)[num.covs]
  covars[,num.covs] <- scale(covars[,num.covs])  # put covariates to mean=0, var=1
  if ( force.linear ) {
    rhs2 <- paste(sapply(num.covs, function(u) { paste('I(', u, '^2)', sep='')}), collapse=' + ')
    fla <- formula(paste(lhs, '~', rhs1, '+', rhs2))
  } else {
    fla <- formula(paste(lhs, rhs1, sep=' ~ '))
  } 
  allowed.fx <- function(degree, pred, parents, namesx, first) {
    bad.idx <- which(namesx %in% no.cross)
    degree == 1 || ! (pred == bad.idx || parents[bad.idx])
  }
  if ( force.linear ) {
    lins <- TRUE
  } else {
    lins <- FALSE
  }
  my.data <- cbind(covars, scale(t(gene.expr[gene.idx,])))  # set genes to mean-0, var-1
  n.terms <- NULL
  eres <- earth(fla, data=my.data, trace=2, degree=3, nk=n.predictors, pmethod='backward', 
                minspan=min.span, endspan=end.span, fast.k=10*n.predictors, linpreds=lins, 
                allowed=allowed.fx, nprune=nprune)
  while ( is.null(n.terms) || n.terms > max.terms ) {
    eres <- update(eres, nprune=nprune)
    coefs <- as.data.frame(sort(apply(eres$coefficients, 1, function(x) { quantile(abs(x), 0.8)}), decreasing=T))
    colnames(coefs) <- 'Beta80Pct'
    if ( 'mean.expr' %in% rownames(coefs) ) {
      drop.idx <- which(rownames(coefs) == 'mean.expr')
      coefs <- coefs[1:(drop.idx-1),,drop=F]
    }
    num.vars <- nrow(coefs)
    is.in.term <- sapply(colnames(covars), function(cname) { grepl(cname, rownames(coefs))})
    terms <- lapply(1:nrow(coefs), function(u) { NULL })
    for ( c.idx in 1:ncol(covars) ) {
      for ( term.idx in which(is.in.term[,c.idx]) ) {
        terms[[term.idx]] <- c(terms[[term.idx]], colnames(covars)[c.idx])
      }
    }
    terms <- lapply(terms, function(s) { paste(s, collapse=':')})
    terms <- unique(unlist(terms))
    if ( '' %in% terms ) {
      terms[terms == ''] <- '1'
    }
    new.fla <- paste('~',paste(sort(terms), collapse=' + '))
    mm <- model.matrix(formula(new.fla), covars)
    n.terms <- ncol(mm)
    nprune = nprune - 1
  }
  list('earth'=eres, 'formula'=new.fla, 'terms'=terms, 'model.matrix'=mm, 'term.imp'=coefs)
}

multiEarthSelect <- function(gene.expr, covars, n.replicates=10, n.cores=4) {
  # run earthSelect `n.replicates` times on 1,000 random samples of genes
  # and return the model importances
  i.seed <- sample.int(10000)
  mclapply(1:n.replicates, function(r) {
    set.seed(i.seed + r)
    earthSelect(gene.expr, covars, n.predictors=200, num.genes=1000, nprune=30, max.terms=20)$term.imp
  }, mc.cores=n.cores)
}
