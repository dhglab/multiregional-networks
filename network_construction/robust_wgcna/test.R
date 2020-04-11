# simple test script for SplicingCoexpression
options(stringsAsFactors=F)
source('SplicingCoexpressionMatrix.R')
test.file <- '/u/nobackup/dhg/chartl/20160725_network_grant/data/GTEx_gene_expression/20170516_V5/20170517.gtex_expression.isoform.brain.good_genes.outlier_rm.txt'
tx.map <- read.table('20171003.transcript_map.txt', header=T)
tx.expr <- read.table(test.file, nrows=25000, header=T)

tx.map$gene <- tx.map$ensembl_gene_id
tx.map$isoform <- tx.map$ensembl_transcript_id

rownames(tx.expr) <- sapply(rownames(tx.expr), function(x) { strsplit(x, '\\.')[[1]][1] })

map.sub <- subset(tx.map, isoform %in% rownames(tx.expr))
expr.sub <- tx.expr[rownames(tx.expr) %in% tx.map$isoform,]

tx.expr.gene1 <-expr.sub[subset(map.sub, gene == 'ENSG00000109670')$isoform,]

tx.expr.gene1.sp <- mat.spearman.transform(as.matrix(tx.expr.gene1), center=T)

#tx.bc <- bicor.transform(as.matrix(expr.sub), map.sub)
#tx.sp <- spearman.transform(as.matrix(expr.sub), map.sub, threads=4)

tx.rmmean <- t(scale(t(as.matrix(expr.sub)), scale=F))

tx.expr.gene1 <- as.matrix(tx.rmmean)[subset(map.sub, gene == 'ENSG00000284448')$isoform,,drop=F]
tx.expr.gene2 <- as.matrix(tx.rmmean)[subset(map.sub, gene == 'ENSG00000012983')$isoform,,drop=F]

tx.expr.gene1 <- semi.orth(tx.expr.gene1)
tx.expr.gene2 <- semi.orth(tx.expr.gene2)

print(tx.expr.gene1[,1:8])
print(tx.expr.gene2[,1:8])

print(tx.expr.gene1 %*% t(tx.expr.gene2))

print(softAdj(svd(tx.expr.gene1 %*% t(tx.expr.gene2))$d))

expr.mat <- as.matrix(tx.rmmean)

adj.mat <- isoformBasedAdjacency(expr.mat, map.sub, ncores=15)

print(adj.mat[1:8,1:8])
print('saving...')
save(adj.mat, 'test25000.Rdata')
