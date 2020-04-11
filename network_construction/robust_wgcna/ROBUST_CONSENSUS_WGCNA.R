library('WGCNA')
code.dir <- Sys.getenv('INTEGRATED_PIPELINE_DIR', '.')
source(paste(code.dir, 'integrated_pipeline_utils.R', sep='/'))

# needed: user.data, TOM.files, min.module.size, cut.height, me.corr, output.dir


args <- commandArgs(trailingOnly = T)
if ( is.null(args) || length(args) != 6 ) {
  stop('USAGE: ROBUST_CONSENSUS_WGCNA.R [datafile] [TOM.list.file] [min.module.size] [cut.height] [module.corr] [output.dir]')
}

#quant <- function(x) { quantile(x, 0.95, type=8) }
quant <- median

user.data <- args[1]
TOM.listfile <- args[2]
min.module.size <- as.integer(args[3])
cut.height <- as.numeric(args[4])
me.corr <- as.numeric(args[5])
output.dir <- args[6]


# marshall output files
consensus.TOM.out = paste(output.dir, 'robust.consensus.TOM.txt', sep='/')
dendro.plot = paste(output.dir, 'robust.consensus.modules.png', sep='/')
rdata.save = paste(output.dir, 'robust.wgcna.Rdata', sep='/')

# have to do this a bit of the consensus at a time
if ( ! file.exists(consensus.TOM.out) ) {
  expected.size = ncol(read.table(user.data, header=T))
  gc()
  expected.size = expected.size*(expected.size+1)/2
  consensus.TOM = rep(0, expected.size)
  print('merging TOM files')
  TOM.files <- read.table(TOM.listfile, header=F, stringsAsFactors=F)[,1]
  CHUNK_SIZE = 1e6
  tom.iters <- lapply(TOM.files, function(tomfile) { fstream(tomfile, double(), header=F, nlines=CHUNK_SIZE)} )
  # read the TOMs into a single matrix
  n.toms <- length(TOM.files)
  consensus.TOM <- rep(0, expected.size)
  chunk = 0
  while ( hasNext(tom.iters[[1]]) ) {
    chunk = 1 + chunk
    idx.end <- min(expected.size, CHUNK_SIZE*chunk)
    idx.start <- 1 + CHUNK_SIZE*(chunk-1)
    tom.tensor <- matrix(0, nrow=1+idx.end-idx.start, ncol=n.toms)
    for ( idx in 1:n.toms ) {
      tom.tensor[,idx] = nextElem(tom.iters[[idx]])
    }
    chunk.consensus <- apply(tom.tensor, 1, quant)
    if ( 1 + idx.end - idx.start != length(chunk.consensus) ) {
      stop('Invalid size!')
    }
    print(paste('start:', idx.start, 'end:', idx.end, 'size:', length(chunk.consensus), 'tom:', length(consensus.TOM)))
    consensus.TOM[idx.start:idx.end] <- chunk.consensus
    print(sprintf('chunk %d (%d/%d = %.3f%%)', chunk, idx.end, expected.size, 100*idx.end/expected.size))
  }
  write.table(consensus.TOM, file=consensus.TOM.out, quote=F, row.names=F, col.names=F)
} else {
  print('reading TOM consensus')
  consensus.TOM <- read.table(consensus.TOM.out)[,1]
  print(length(consensus.TOM))
}
consensus.TOM <- as.matrix(wrap.symm(consensus.TOM))
tree = hclust(as.dist(1-consensus.TOM), method='average')
modules = cutreeHybrid(dendro=tree, distM=1-consensus.TOM, deepSplit=2, pamStage=F,
                        pamRespectsDendro=F, minClusterSize=min.module.size, cut.height)
modules$colors <- labels2colors(modules$labels)
print(table(modules$colors))
data <- read.table(user.data, header=T)
print(length(modules$colors))
print(dim(data))
modules.merged <- mergeCloseModules(data, modules$labels, cutHeight=1-me.corr)
print(names(modules.merged))

rownames(consensus.TOM) <- colnames(data)

save(tree, modules, modules.merged, file=rdata.save)
save(consensus.TOM, file=paste(rdata.save, 'named.TOM.Rdata', sep='.'))

png(dendro.plot, height=800, width=1200)
plotDendroAndColors(tree, cbind('merged'=modules.merged$colors, 'orig'=modules$colors), addGuide=T, dendroLabels=F)
dev.off()

cuts <- 1 - exp(seq(log(0.2), log(0.001), length.out=20))
idx <- 1
for ( cut in cuts ) {
  mods <- cutreeHybrid(dendro=tree, distM=1-consensus.TOM, cut, minClusterSize=min.module.size)
  mods.merged <- mergeCloseModules(data, mods$labels)
  m.kme <- mods.merged$newMEs
  m.colors <- labels2colors(mods.merged$colors)
  if ( length(unique(mods$labels)) < 5 ) {
    next
  }
  mods.merged <- tryCatch({mergeCloseModules(data, mods$labels)}, error=function(cond){ NA })
  if ( is.na(mods.merged) ) {
    next
  }
  m.kme <- mods.merged$newMEs
  m.colors <- labels2colors(mods.merged$colors)
  if ( length(unique(m.colors)) < 3 ) {
    next
  }
  m.colors.orig <- labels2colors(mods$labels)
  png(sprintf('modules.cut.%d.png', idx), height=800, width=1200)
  plotDendroAndColors(tree, cbind('merged'=m.colors, 'orig'=m.colors.orig), addGuide=T, dendroLabels=F, main=cut)
  dev.off()
  idx <- 1 + idx
}

