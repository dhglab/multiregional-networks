library(Rtsne)
library(pdist)
options(stringsAsFactors=F)
load('module_eigengenes.withlandmarks.Rda')
for ( tissue in tissues ) {
  print(sprintf('Plotting %s', tissue))
  expr <- read.table(expression.data[tissue], header=T, sep='\t')
  expr.lnd <- t(expr[landmark.genes,])
  mods.to.use <- c(tissue, parents.manual[[tissue]])
  print(mods.to.use)
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

