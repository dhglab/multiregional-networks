# GLOBAL MERGING
options(stringsAsFactors = F)
library(abind)
library(fastcluster)
rm(list=ls())
load('./module_eigengenes.withlandmarks.Rda')
module.defs <- read.table('module.table.txt', header=T, row.names=1)
module.color.defs <- read.table('module.color.table.txt', header=T, row.names=1)
module.defs.orig <- module.defs
for ( col in 1:ncol(module.defs) ) {
  module.defs[,col] <- sapply(module.defs[,col], function(q) { 
     r <- strsplit(q, '\\.')[[1]]
     paste('.', r[length(r)], sep='')
  })
}
module.inds <- model.matrix(~ . - 1, data=module.defs)
overlap.counts <- t(module.inds) %*% module.inds
module.ind.jac <- matrix(0, nrow=ncol(module.inds), ncol=ncol(module.inds))
for ( idx in 1:(ncol(module.inds)-1) ) {
  print(idx)
  for ( idx2 in (1+idx):ncol(module.inds) ) {
    module.ind.jac[idx, idx2] = overlap.counts[idx, idx2]/min(overlap.counts[idx, idx], overlap.counts[idx2, idx2]) #(overlap.counts[idx,idx] + overlap.counts[idx2, idx2] - overlap.counts[idx, idx2])
    module.ind.jac[idx2, idx] <- module.ind.jac[idx, idx2]
  }
}
diag(module.ind.jac) <- 1
jaccard.dissim <- 1 - module.ind.jac
rownames(jaccard.dissim) <- colnames(jaccard.dissim) <- sprintf('ME%s',colnames(module.inds))
M0s <- which(grepl('M0', colnames(jaccard.dissim)))
jaccard.dissim <- jaccard.dissim[-M0s,-M0s]

Q <- abind(eigen.corrs, along=3)
consensus.dissim <- apply(1 - Q, c(1,2), function(x) { quantile(x, 0.99, na.rm=T)})
colnames(consensus.dissim) <- sapply(colnames(consensus.dissim), function(x) {
  if ( grepl('NS.SCTX', x) ) {
    q <- strsplit(x, '\\.')[[1]]
    sprintf('MENS.SCTX.%s', q[length(q)])
  } else {
    q <- strsplit(x, '\\.')[[1]]
    sprintf('%s.%s', q[1], q[3])
  }
})  # MEBRNACC.BRNACC.M1
rownames(consensus.dissim) <- colnames(consensus.dissim)
M0s <- which(grepl('M0', colnames(consensus.dissim)))
consensus.dissim <- consensus.dissim[-M0s,-M0s]
print(head(colnames(jaccard.dissim)))
print(head(colnames(consensus.dissim)))
jaccard.dissim <- jaccard.dissim[colnames(consensus.dissim), colnames(consensus.dissim)]
consensus.dissim[is.na(consensus.dissim)] <- 1 # jaccard.dissim[is.na(consensus.dissim)] 
consensus.dissim.wjac <- (2 * consensus.dissim + 6 * jaccard.dissim)/8
consensus.dissim.orig <- consensus.dissim.wjac

orig.names <- colnames(consensus.dissim.orig)
merged.names <- colnames(consensus.dissim.orig)
merged.mods <- NULL
any.merged <- T
ii <- 1
CUT_HEIGHT=0.35
while ( any.merged ) {
  print(sprintf('merge %d', ii))
  ii <- 1 + ii
  any.merged <- F
  if ( is.null(merged.mods) ) {
    merged.names <- unique(merged.names)
  } else {
    merged.names <- unique(merged.mods$merged.name)
  }
  consensus.dissim <- consensus.dissim.orig[merged.names, merged.names]
  orig.names <- merged.names
  tree <- hclust(as.dist(consensus.dissim), 'average')
  plot(tree)
  branches <- as.factor(WGCNA::moduleNumber(dendro=tree, cutHeight=CUT_HEIGHT, minSize=1))
  uniqueBranches <- levels(branches)
  nBranches <- nlevels(branches)
  numOnBranch <- table(branches)

  for ( branch in 1:nBranches ) {  
    if ( numOnBranch[branch] > 1 ) {
      modsOnThisBranch <- names(branches)[branches == uniqueBranches[branch]]
      for ( j in 2:length(modsOnThisBranch) ) {
        any.merged <- T
        merged.names[merged.names == modsOnThisBranch[j]] <- modsOnThisBranch[1]
      }
    }
  }
  if ( is.null(merged.mods) ) {
    merged.mods <- data.frame(orig.name=orig.names, merged.name=merged.names)
  } else {
    next.merge <- data.frame(merged.name=orig.names, merged.final=merged.names)
    mmods <- merge(merged.mods, next.merge, on='merged.name')
    merged.mods <- data.frame(orig.name=mmods$orig.name, merged.name=mmods$merged.final)
  }
  td.names <- c('WHOLE_BRAIN', 'ALL', 'CEREBELLUM', 'MECTX', 'SUBCTX', 'NS\\.SCTX', 'STR', 'BGA', 'BROD')
  for ( mn in unique(merged.mods$merged.name) ) {
    mm.sub <- subset(merged.mods, merged.name == mn)
    td.idx <- which(sapply(td.names, function(x) { grepl(x, mn)}))
    if ( length(td.idx) == 0 ) {
      td.idx <- length(td.names)
    } else {
      td.idx <- td.idx - 1
    }
    if ( td.idx == 0 ) {
      next
    }
    # take the top-level name if possible
    for ( td.name in td.names[1:td.idx]) {
      if ( any(grepl(td.name, mm.sub$orig.name)) ) {
        newname <- mm.sub$orig.name[which(grepl(td.name, mm.sub$orig.name))[1]]
        merged.mods$merged.name[merged.mods$merged.name == mn] = newname
        break
      }
    }
  }
}
print('Done merging! fixing order')
print(table(merged.mods$merged.name))
# now go back to the original module definitions to define final names
module.lists <- list()
mnames.fixed.odr <- unique(merged.mods$merged.name)
for ( modname in mnames.fixed.odr ) {
  super.tissue <- gsub('^ME', '', strsplit(modname, '\\.M')[[1]][1])
  mod.num <- strsplit(modname, '\\.M')[[1]]
  mod.num <- mod.num[length(mod.num)]
  mod.sub <- subset(merged.mods, merged.name == modname)
  tissue.modname <- sprintf('%s.M%s', super.tissue, mod.num)
  genes <- rownames(module.defs.orig)[module.defs.orig[,super.tissue] == tissue.modname]
  module.lists[[super.tissue]][[modname]] <- genes
}

print('step1 done')

for ( tis in names(module.lists) ) {
  mods <- module.lists[[tis]]
  sizes <- sapply(mods, length)
  names <- sprintf("M%d", 1:length(sizes))
  orig.names <- names(mods)
  for ( idx in 1:length(orig.names) ) {
    name <- orig.names[idx]
    newname <- sprintf('%s.%s', tis, names[idx])
    merged.mods$merged.name[merged.mods$merged.name == name] <- newname
  }
}

print('step2 done')

write.table(merged.mods, file='consensus_merged_modules.txt', quote=F)

tissues <- sapply(merged.mods$orig.name, function(x) { gsub('^ME', '', strsplit(x, '\\.M')[[1]][1]) })
mod.nums <- sapply(merged.mods$orig.name, function(x) {
  mn <- strsplit(x, '\\.M')[[1]]
  as.integer(mn[length(mn)])
})

save(list=ls(), file='debug.Rdata')
module.merging <- data.frame(tissue=tissues, module.label=sprintf("M%d", mod.nums), final.module=merged.mods$merged.name)
rownames(module.merging) <- NULL
write.table(module.merging, file='consensus_merged_modules.mapping.txt')
dev.off()
