library(Rtsne)
library(pdist)
load('module_eigengenes.withlandmarks.Rda')
options(stringsAsFactors=F)
me.merged <- lapply(1:length(module.eigengenes), function(idx) {
  M <- do.call(cbind, module.eigengenes[[idx]])
  rownames(M) <- sprintf('%s.%d', names(module.eigengenes)[idx], 1:nrow(M))
  M
})
me.merged <- do.call(rbind, me.merged)

pdf('cross_tissue.module.merging.pdf')
mods.to.use <- c(tissues, names(tissue.hierarchy))
colnames.modules <- sapply(colnames(me.merged), function(s) {
    if ( grepl('NS.SCTX', s) ) {
      'NS.SCTX'
    } else {
      strsplit(s, '\\.')[[1]][1]
    }
  })
mod.colors <- rainbow(length(mods.to.use), alpha=0.8)[factor(colnames.modules, levels=mods.to.use)]
print(head(mod.colors))
tsne.res2 <- Rtsne(t(me.merged), theta=0.0, max_iter=1000, perplexity=8)
rownames(tsne.res2$Y) <- colnames(me.merged)
plot(tsne.res2$Y, col=mod.colors, pch=16)
legend('topleft', legend=mods.to.use, col=rainbow(length(mods.to.use), alpha=0.8)[factor(mods.to.use, levels=mods.to.use)], pch=16)
plot(tsne.res2$Y, col=mod.colors, pch=16)
text(tsne.res2$Y, colnames(me.merged), cex=0.2)


# build the top/down and bottom/up hierarchies
edge.df <- data.frame()
for ( tissue in tissues ) {
  td.hier <- c(rev(parents.manual[[tissue]]), tissue)  # whole brain -> all -> subctx -> ns.sctx -> brnacc
  for ( idx in 1:(length(td.hier)-1) ) {
    print(td.hier[idx])
    print(td.hier[1+idx])
    p.idx <- grepl(sprintf('^%s',td.hier[idx]), colnames(me.merged))
    c.idx <- grepl(sprintf('^%s',td.hier[1+idx]), colnames(me.merged))
    r.idx <- grepl(sprintf('^%s', tissue), rownames(me.merged))
    XX = scale(me.merged[r.idx,p.idx])
    YY = scale(me.merged[r.idx,c.idx])
    CC <- WGCNA::bicor(XX,YY)
    D <- as.matrix(pdist(X=t(XX), Y=t(YY)))
    D2 <- as.matrix(pdist(tsne.res2$Y[colnames(XX),], tsne.res2$Y[colnames(YY),]))
    D2 <- D
    rownames(D2) <- colnames(XX)
    colnames(D2) <- colnames(YY)
    heatmap(CC, cexRow=0.4, cexCol=0.4)
    td.conns <- apply(D2, 1, which.min)
    bu.conns <- apply(D2, 2, which.min)
    for ( i in 1:nrow(D2) ) {
      for ( j in 1:ncol(D2) ) {
        if ( CC[i, j] > 0.8 ) {
          entry <- data.frame(mod.source=rownames(CC)[i], mod.destination=colnames(CC)[j], edge.type='cor', distance=D2[i, j], dist.source=tissue, cor=CC[i,j])
          edge.df <- rbind(edge.df, entry)
        }
      }
      entry <- data.frame(mod.source=rownames(D2)[i], mod.destination=colnames(D2)[td.conns[i]], edge.type='TD', distance=D2[i, td.conns[i]], dist.source=tissue, cor=CC[i, td.conns[i]])
      edge.df <- rbind(edge.df, entry)
    }
    for ( j in 1:ncol(D2) ) {
      entry <- data.frame(mod.source=colnames(D2)[j], mod.destination=rownames(D2)[bu.conns[j]], edge.type='BU', distance=D2[bu.conns[j], j], dist.source=tissue, cor=CC[bu.conns[j], j])
      edge.df <- rbind(edge.df, entry)
    }
  }
}


write.table(edge.df, quote=F, file='all_tissue.closest_modules.txt')

plot(tsne.res2$Y, col=mod.colors, pch=16)
for ( idx in 1:nrow(edge.df) ) {
    print(edge.df[idx,])
    pt.start <- tsne.res2$Y[edge.df$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df$mod.destination[idx],]
    arr.col <- if ( edge.df$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
}
  edge.df2 <- edge.df
  edge.df2 <- subset(edge.df2, cor > 0.8)
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
  write.table(edge.df2, quote=F, file='alltissue.closest_modules.outlier_edge_rm.txt')
  plot(tsne.res2$Y, col=mod.colors, pch=16)
  for ( idx in 1:nrow(edge.df2) ) {
    print(edge.df2[idx,])
    pt.start <- tsne.res2$Y[edge.df2$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df2$mod.destination[idx],]
    arr.col <- if ( edge.df2$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
  }

  plot(tsne.res2$Y, col=mod.colors, pch=16)
  for ( idx in 1:nrow(edge.df2) ) {
    print(edge.df2[idx,])
    pt.start <- tsne.res2$Y[edge.df2$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df2$mod.destination[idx],]
    if ( edge.df2$edge.type[idx] == 'BU' ) { next }
    arr.col <- if ( edge.df2$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
  }

  plot(tsne.res2$Y, col=mod.colors, pch=16)
  for ( idx in 1:nrow(edge.df2) ) {
    print(edge.df2[idx,])
    pt.start <- tsne.res2$Y[edge.df2$mod.source[idx],]
    pt.end <- tsne.res2$Y[edge.df2$mod.destination[idx],]
    if ( edge.df2$edge.type[idx] == 'TD' ) { next }
    arr.col <- if ( edge.df2$edge.type[idx] == 'TD' ) '#555555AA' else '#FF00FFAA'
    arrows(x0=pt.start[1], x1=pt.end[1], y0=pt.start[2], y1=pt.end[2], code=2, length=0.05, col=arr.col)
  }
dev.off()


                    
