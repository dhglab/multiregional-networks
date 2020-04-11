options(stringsAsFactors=F)
source('SplicingCoexpressionMatrix.R')
test.file <- '/u/nobackup/dhg/chartl/20160725_network_grant/data/GTEx_gene_expression/20170516_V5/20170517.gtex_expression.isoform.brain.good_genes.outlier_rm.txt'
tx.map <- read.table('20171003.transcript_map.txt', header=T)
tx.expr <- read.table(test.file, header=T)

tx.map$gene <- tx.map$ensembl_gene_id
tx.map$isoform <- tx.map$ensembl_transcript_id

rownames(tx.expr) <- sapply(rownames(tx.expr), function(x) { strsplit(x, '\\.')[[1]][1] })

map.sub <- subset(tx.map, isoform %in% rownames(tx.expr))
expr.sub <- tx.expr[rownames(tx.expr) %in% tx.map$isoform,]

tissues <- sapply(colnames(tx.expr), function(x) { strsplit(x, '\\.')[[1]][3] })

print(head(tissues))

tx.expr.sub <- tx.expr[,tissues == 'BRNCTX']

print(dim(tx.expr.sub))

#tx.bc <- bicor.transform(as.matrix(expr.sub), map.sub)
#tx.sp <- spearman.transform(as.matrix(expr.sub), map.sub, threads=4)

tx.rmmean <- t(scale(t(as.matrix(expr.sub)), scale=F))

expr.mat <- as.matrix(tx.rmmean)

adj.mat <- isoformBasedAdjacency(expr.mat, map.sub, ncores=1, pow=10, njobs=900)

save(adj.mat, file='20171010_isoform_gene_adj.BRNCTX.Rdata')

