library(Rtsne)
library(WGCNA)
options(stringsAsFactors=F)
tissues <- c('BRNACC', 'BRNAMY', 'BRNCBH', 'BRNCBL', 'BRNCDT', 'BRNCTX', 'BRNCTXB24', 'BRNCTXBA9', 'BRNHIP', 'BRNHYP', 'BRNPUT', 'BRNSNA')
tissue.hierarchy <- list(
               'BROD'=c('BRNCTXB24', 'BRNCTXBA9'), 
               'BGA'=c('BRNCDT', 'BRNPUT'), 
               'STR'=c('BRNACC', 'BGA'),
               'NS.SCTX'=c('BRNAMY', 'BRNHIP', 'BRNHYP', 'BRNSNA'), 
               'SUBCTX'=c('STR', 'NS.SCTX'),
               'CTX'=c('BRNCTX', 'BROD'), 'ALL'=c('CTX', 'SUBCTX'), 
               'CEREBELLUM'=c('BRNCBL', 'BRNCBH'),
               'WHOLE_BRAIN'=c('ALL','CEREBELLUM'))

td.recurse <- function(x) {
  if ( x %in% names(tissue.hierarchy) ) {
    c(x, do.call(c, lapply(tissue.hierarchy[[x]], td.recurse)))
  } else {
    x
  }
}

parents.manual <- list(
 'BRNACC' = c('STR', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNAMY' = c('NS.SCTX', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNCBH' = c('CEREBELLUM', 'WHOLE_BRAIN'),
 'BRNCBL' = c('CEREBELLUM', 'WHOLE_BRAIN'),
 'BRNCDT' = c('BGA', 'STR', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNCTX' = c('CTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNCTXB24' = c('BROD', 'CTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNCTXBA9' = c('BROD', 'CTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNHIP' = c('NS.SCTX', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNHYP' = c('NS.SCTX', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNPUT' = c('BGA', 'STR', 'SUBCTX', 'ALL', 'WHOLE_BRAIN'),
 'BRNSNA' = c('NS.SCTX', 'SUBCTX', 'ALL', 'WHOLE_BRAIN')
)

expression.data <- sapply(tissues, function(q) { sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/%s/20171010.gtex_expression.%s.lm_regressed.txt', q, q)})
print(expression.data)

module.data <- sapply(tissues, function(q) { sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/%s/rWGCNA/robust.wgcna.bs.cut.Rdata', q) })
module.data <- c(module.data, sapply(names(tissue.hierarchy), function(q) { sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/hierarchical/%s.min_cor.cut.Rdata', q)}))
tom.data <- sapply(tissues, function(q) { sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/%s/rWGCNA/robust.wgcna.Rdata.named.TOM.Rdata', q) })
tom.data <- c(tom.data, sapply(names(tissue.hierarchy), function(q) { sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/hierarchical/%s.consensus.TOM.Rda', q)}))

modules.by.tom <- lapply(names(module.data), function(dname) {
  dfile <- module.data[[dname]]
  env <- new.env()
  load(dfile, env=env)
  print(dname)
  if ( 'module.labels' %in% names(env$modules) ) {
    x <- env$modules$module.labels
    x[is.na(x)] <- 'M0'
    x <- sprintf('%s.%s', dname, x)
    names(x) <- names(env$modules$module.colors)
  } else {
    x <- sprintf('%s.M%d', dname, env$modules$labels)
  }
  print(head(x))
  if ( any(is.na(x)) ) {
    stop(sprintf('NA values in %s (%s)', dname, dfile))
  }
  sprintf('%s.%s', dname, x)
})
names(modules.by.tom) <- names(module.data)

colors.by.tom <- lapply(names(module.data), function(dname) {
  dfile <- module.data[[dname]]
  env <- new.env()
  load(dfile, env=env)
  print(dname)
  if ( 'module.colors' %in% names(env$modules) ) {
    x <- env$modules$module.colors
    x[is.na(x)] <- 'grey'
    names(x) <- names(env$modules$module.colors)
  } else {
    x <- env$modules$module.colors
  }
  x
})
names(colors.by.tom) <- names(module.data)

for ( tissue in names(colors.by.tom) ) {
  png(sprintf('20171204.%s.tree.png', tissue), height=800, width=1400)
  t.env <- new.env()
  load(tom.data[[tissue]], env=t.env)
  tom <- t.env[[names(t.env)[1]]]
  rm(t.env)
  tree <- hclust(as.dist(1-tom), 'average')
  dend.colors <- list()
  if ( tissue %in% tissues ) {
    # bottom-up colors
    dend.colors[[tissue]] <- colors.by.tom[[tissue]]
    for ( pa.tis in parents.manual[[tissue]] ) {
      dend.colors[[pa.tis]] <- colors.by.tom[[pa.tis]]
    }
  } else {
    # top-down colors
    dend.colors[[tissue]] <- colors.by.tom[[tissue]]
    for ( ch.tis in td.recurse(tissue) ) {
      dend.colors[[ch.tis]] <- colors.by.tom[[ch.tis]]
    }
  }
  cnames <- names(dend.colors)
  dend.colors <- do.call(cbind, dend.colors)
  plotDendroAndColors(tree, dend.colors, groupLabels=cnames, dendroLabels=F)
  dev.off()
}

all.modules <- do.call(cbind, modules.by.tom)
all.colors <- do.call(cbind, colors.by.tom)
rownames(all.modules) <- rownames(read.table(expression.data[1], header=T))
rownames(all.colors) <- rownames(read.table(expression.data[1], header=T))
write.table(all.modules, file='module.table.txt', quote=F)
write.table(all.colors, file='module.color.table.txt', quote=F)

getLandmarks <- function(expr, n.landmarks) {
  landmark.idx <- sample.int(nrow(expr), 1)
  u <- t(t(expr) - as.numeric(expr[landmark.idx,]))
  dist.mat <- matrix(rowSums(u*u), ncol=1)
  for ( nxt in 1:(n.landmarks - 1) ) {
    #print(sprintf('%d of %d', 1 + nxt, n.landmarks))
    md <- apply(dist.mat, 1, min)
    next.idx <- which.max(md)
    u <- t(t(expr) - as.numeric(expr[next.idx,]))
    dist.mat <- cbind(md, rowSums(u*u))  # only care about minimum distance
    landmark.idx <- c(landmark.idx, next.idx)
  }
  landmark.idx
}


n.tissue.landmarks <- 65
landmark.genes <- NULL
module.eigengenes <- list()

for ( tissue in tissues ) {
  expr <- read.table(expression.data[tissue], header=T, sep='\t')
  expr.scl <- t(scale(t(expr)))
  landmark.genes <- c(landmark.genes, rownames(expr)[getLandmarks(expr.scl, n.tissue.landmarks)])
  module.eigengenes[[tissue]] <- lapply(names(modules.by.tom), function(tname) {
    eigens <- WGCNA::moduleEigengenes(t(expr), modules.by.tom[[tname]], grey=sprintf('%s.M0', tname), verbose=2)
    print(eigens$eigengenes[1:4,1:4])
    eigens$eigengenes
  })
  names(module.eigengenes[[tissue]]) <- names(modules.by.tom)
}

landmark.genes <- unique(landmark.genes)

library(pdist)
eigen.corrs <- list()
for ( tissue in tissues ) {
  print(sprintf('Plotting %s', tissue))
  expr <- read.table(expression.data[tissue], header=T, sep='\t')
  expr.lnd <- t(expr[landmark.genes,])
  mods.to.use <- c(tissue, parents.manual[[tissue]])
  print(mods.to.use)
  print(names(module.eigengenes[[tissue]]))
  eigen.data <- do.call(cbind, module.eigengenes[[tissue]][mods.to.use])
  print(eigen.data[1:4,1:4])
  vectors <- t(scale(cbind(expr.lnd, eigen.data)))
  plot.colors <- rep('#AAAAAA44', ncol(expr.lnd))
  colnames.modules <- sapply(colnames(eigen.data), function(s) { 
    if ( grepl('NS.SCTX', s) ) {
      'NS.SCTX'
    } else {
      strsplit(s, '\\.')[[1]][1]
    }
  })
  print(mods.to.use)
  print(head(colnames.modules))
  mod.colors <- rainbow(length(mods.to.use), alpha=0.8)[factor(colnames.modules, levels=mods.to.use)]
  print(head(mod.colors))
  plot.colors <- c(plot.colors, mod.colors)
  tsne.res <- Rtsne(vectors, theta=0.0, max_iter=1000, perplexity=8)
  pdf(sprintf('%s.hierarchy.tsne.pdf', tissue))
  plot(tsne.res$Y, col=plot.colors, pch=16)
  legend('topleft', legend=mods.to.use, col=rainbow(length(mods.to.use), alpha=0.8)[factor(mods.to.use)], pch=16)
  mod.labels <- c(rep('', ncol(expr.lnd)), colnames(eigen.data))
  plot(tsne.res$Y, col=plot.colors, pch=16)
  text(tsne.res$Y, mod.labels, cex=0.2)
  tsne.res2 <- Rtsne(t(eigen.data), theta=0.0, max_iter=1000, perplexity=8)
  rownames(tsne.res2$Y) <- sapply(colnames(eigen.data), function(x) { 
   if ( grepl('NS.SCTX', x) ) {
     paste(strsplit(x, '\\.')[[1]][-c(1,2)], collapse='.')
   } else {
     paste(strsplit(x, '\\.')[[1]][-1], collapse='.')
   }
  })
  plot(tsne.res2$Y, col=mod.colors, pch=16)
  legend('topleft', legend=mods.to.use, col=rainbow(length(mods.to.use), alpha=0.8)[factor(mods.to.use, levels=mods.to.use)], pch=16)
  plot(tsne.res2$Y, col=mod.colors, pch=16)
  text(tsne.res2$Y, mod.labels, cex=0.2)

  # build the top/down and bottom/up hierarchies
  td.hier <- rev(mods.to.use)  # whole brain -> all -> subctx -> ns.sctx -> brnacc
  edge.df <- data.frame()
  na.mes <- lapply(names(module.eigengenes[[tissue]]), function(pname) {
    eg <- module.eigengenes[[tissue]][[pname]]
    if ( pname %in% mods.to.use ) {
      eg
    } else {
      rn <- rownames(eg)
      cn <- colnames(eg)
      eg <- matrix(NA, nrow=nrow(eg), ncol=ncol(eg))
      rownames(eg) <- rn
      colnames(eg) <- cn
      eg
    }
  })
  eigen.corrs[[tissue]] <- bicor(do.call(cbind, na.mes), use='pairwise.complete.obs')
  for ( idx in 1:(length(td.hier)-1) ) {
    print(td.hier[idx])
    print(td.hier[1+idx])
    XX = scale(module.eigengenes[[tissue]][[td.hier[idx]]])
    YY = scale(module.eigengenes[[tissue]][[td.hier[idx+1]]])
    CC <- cor(XX,YY)
    D <- as.matrix(pdist(X=t(XX), Y=t(YY)))
    D2 <- as.matrix(pdist(tsne.res2$Y[colnames(XX),], tsne.res2$Y[colnames(YY),]))
    rownames(D2) <- colnames(XX)
    colnames(D2) <- colnames(YY)
    heatmap(CC, cexRow=0.4, cexCol=0.4)
    td.conns <- apply(D2, 1, which.min)
    bu.conns <- apply(D2, 2, which.min)
    for ( i in 1:nrow(D2) ) {
      entry <- data.frame(mod.source=rownames(D2)[i], mod.destination=colnames(D2)[td.conns[i]], edge.type='TD', distance=D2[i, td.conns[i]])
      edge.df <- rbind(edge.df, entry)
    }
    for ( j in 1:ncol(D2) ) {
      entry <- data.frame(mod.source=colnames(D2)[j], mod.destination=rownames(D2)[bu.conns[j]], edge.type='BU', distance=D2[bu.conns[j], j])
      edge.df <- rbind(edge.df, entry)
    }
  }
  write.table(edge.df, quote=F, file=sprintf('%s.closest_modules.txt', tissue))

  plot(tsne.res2$Y, col=mod.colors, pch=16)
  for ( idx in 1:nrow(edge.df) ) {
    print(edge.df[idx,])
    pt.start <- tsne.res2$Y[edge.df$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df$mod.destination[idx],]
    arr.col <- if ( edge.df$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
  }
  edge.df2 <- edge.df
  edge.outlier <- T
  ZT <- 3.0
  while ( edge.outlier ) {
    hist(edge.df2$distance, 25, col='cornflowerblue')
    points(rep(mean(edge.df2$distance) + ZT*sd(edge.df2$distance), 2), c(-1, 100), type='l', lwd=2, lty=2, col='magenta')
    edge.outlier <- F
    edge.df2$dist.z <- scale(edge.df2$distance)
    edge.outlier <- any(edge.df2$dist.z > ZT)
    edge.df2 <- subset(edge.df2, abs(dist.z) < ZT)
  }
  write.table(edge.df2, quote=F, file=sprintf('%s.closest_modules.outlier_edge_rm.txt', tissue))
  plot(tsne.res2$Y, col=mod.colors, pch=16)
  for ( idx in 1:nrow(edge.df2) ) {
    print(edge.df2[idx,])
    pt.start <- tsne.res2$Y[edge.df2$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df2$mod.destination[idx],]
    arr.col <- if ( edge.df2$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
  }
  dev.off()
}

save(list=ls(), file='module_eigengenes.withlandmarks.Rda')
pdf('all_tissue.consensus.pdf')
ecorr.vec <- do.call(c, lapply(eigen.corrs, function(X) { X[lower.tri(X)] }))
quantile.match <- function(source.vec, target.vec) {
  sv <- sort(source.vec, decreasing=F)
  tar.ranks <- as.integer(rank(target.vec, na.last='keep') * length(sv)/sum(! is.na(target.vec)))
  sv[tar.ranks]
}

print('matching...')
eigen.corrs.match <- lapply(eigen.corrs, function(X) {
  Y <- 0 * X
  ym <- quantile.match(source.vec=X[lower.tri(X)], target.vec=ecorr.vec)
  Y[lower.tri(Y)] <- ym
  Y <- t(Y)
  Y[lower.tri(Y)] <- ym
  diag(Y) <- 1
  Y
})
print('dissim...')
tissue.dissims <- lapply(eigen.corrs.match, function(X) {
  1 - X
})
tissue.dists <- lapply(eigen.corrs.match, acos)
library(abind)
library(WGCNA)
print('consensus...')
consensus.dissim <- apply(abind(tissue.dissims, along=3), c(1,2), max, na.rm=T)
print('rename...')
print(dim(consensus.dissim))
print(dim(eigen.corrs[[1]]))
print(dim(eigen.corrs[[2]]))
rownames(consensus.dissim) <- colnames(consensus.dissim) <- colnames(eigen.corrs[[1]])
consensus.dist <- apply(abind(tissue.dists, along=3), c(1,2), max, na.rm=T)
#consensus.tsne <- Rtsne(consensus.dist, is_distance=TRUE, perplexity=6, theta=0.05)
#plot(consensus.tsne$Y)
print('clustering...')
print(consensus.dissim[1:5,1:5])
print(typeof(consensus.dissim))
write.table(consensus.dissim, file='consensus.dissim.txt')
tree <- hclust(as.dist(consensus.dissim), method='average')
print('done...')
plot(tree)
branches <- as.factor(moduleNumber(dendro=tree, cutHeight=0.25, minSize=1))
uniqueBranches <- levels(branches)
nBranches <- nlevels(branches)
numOnBranch <- table(branches)
merged.names <- colnames(consensus.dissim)
for ( branch in 1:nBranches ) {
  if ( numOnBranch[branch] > 1 ) {
    modsOnThisBranch <- names(branches)[branches == uniqueBranches[branch]]
    for ( j in 2:length(modsOnThisBranch) ) {
      merged.names[merged.names == modsOnThisBranch[j]] <- modsOnThisBranch[1]
    }
  }
}
merged.mods <- data.frame(orig.name=colnames(consensus.dissim), merged.name=merged.names)
write.table(merged.mods, file='consensus_merged_modules.txt', quote=F)
dev.off()

