library(Rtsne)
load('module_eigengenes.withlandmarks.Rda')
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

