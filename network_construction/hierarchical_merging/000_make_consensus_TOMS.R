library('WGCNA')
# brain consensus modules

TISSUES <- c('BRNACC', 'BRNAMY', 'BRNCBH', 'BRNCBL', 'BRNCDT', 'BRNCTX', 'BRNCTXB24', 'BRNCTXBA9', 'BRNHIP', 'BRNHYP', 'BRNPUT', 'BRNSNA')
MERGES <- list('BROD'=c('BRNCTXB24', 'BRNCTXBA9'), 'BGA'=c('BRNCDT', 'BRNPUT'), 'STR'=c('BRNACC', 'BGA'),
               'NS.SCTX'=c('BRNAMY', 'BRNHIP', 'BRNHYP', 'BRNSNA'), 'SUBCTX'=c('STR', 'NS.SCTX'),
               'CTX'=c('BRNCTX', 'BROD'), 'ALL'=c('CTX', 'SUBCTX'), 'CEREBELLUM'=c('BRNCBL', 'BRNCBH'),
               'WHOLE_BRAIN'=c('ALL','CEREBELLUM'))
               


tmp.env <- new.env()
tissue.toms <- list()

t.idx <- 1
for ( tis in TISSUES ) {
  tom.file <- sprintf('/u/nobackup/dhg/chartl/20160725_network_grant/outputs/GTEx/20171024_v5.5_networks/%s/rWGCNA/robust.wgcna.Rdata.named.TOM.Rdata', tis)
  load(tom.file, envir=tmp.env)
  tissue.toms[[t.idx]] <- tmp.env[['consensus.TOM']]
  print(tissue.toms[[t.idx]][1:5,1:5])
  t.idx <- 1 + t.idx
}

numGenes <- nrow(tissue.toms[[1]])
gene.names <- rownames(tissue.toms[[1]])
print(head(gene.names))
tissue.toms <- do.call(cbind, lapply(tissue.toms, as.dist))
colnames(tissue.toms) <- TISSUES

print('tissue.tom shape')
print(dim(tissue.toms))

toms2info <- function(tomlist) {
  print(tomlist)
  print(colnames(tissue.toms))
  list(TOMSimilarities     = list(toms=tissue.toms[,tomlist]),
       blocks              = rep(1, numGenes),
       blockGenes          = list(toms=1:numGenes),
       saveTOMs            = F,
       goodSamplesAndGenes = list(allOK=T, goodGenes=rep(T, numGenes)),
       nGGenes             = numGenes,
       gBlocks             = rep(1, numGenes),
       nSets               = length(tomlist),
       setNames            = colnames(tissue.toms)[tomlist])
}


runConsensus <- function(tis.vec) {
  tominfo <- toms2info(tis.vec)
  print(dim(tominfo$TOMSimilarities[[1]]))
  ct <- consensusTOM(individualTOMInfo   = tominfo, 
                     networkCalibration  = 'full quantile',
		     calibrationQuantile = 0.95,  # default but ignored
                     #consensusQuantile   = 0.2,  # effectively the min
                     saveConsensusTOMs   = F,
                     useMean             = T,
                     returnTOMs          = T, 
                     useDiskCache        = F, 
                     cacheDir            = '/u/scratch/c/chartl/', 
                     nThreads            = 1, 
                     verbose             = T)
  tom <- as.matrix(ct$consensusTOM[[1]])
  rownames(tom) <- gene.names
  colnames(tom) <- gene.names
  tom
}

for ( idx in 1:length(MERGES) ) {
  mg.name <- names(MERGES)[idx]
  print(sprintf('Merging %s', mg.name))
  mg.defs <- MERGES[[idx]]
  cons.tom <- runConsensus(mg.defs)
  save(cons.tom, file=sprintf('%s.consensus.TOM.Rda', mg.name))
  tt.names <- colnames(tissue.toms)
  tissue.toms <- cbind(tissue.toms, as.dist(cons.tom))
  colnames(tissue.toms) <- c(tt.names, mg.name)
}
  
  
