## Script to set up data for ARACNE-AP
## Gokul Ramaswami 10-10-2016

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# data paths
datadir <- "../data/"

# total brain network
datExpr <- read.table(paste0(datadir, "Brain-AllTissues.samples_removed.regressed.ensid.txt"), header=TRUE, row.names=1) # datExpr
TFs <- read.table(paste0(datadir, "Intersect_GO_0003677_0006355_TFs.txt"), header=FALSE) # List of Transcription Factors based on GO
TFs_overlap <- TFs$V1[TFs$V1 %in% rownames(datExpr)]
# write out to files
write.table(as.data.frame(TFs_overlap), file=paste0(datadir, "TFs_totalBrain.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
# fix line 1 with sample names
datExpr_tmp <- rbind(c(colnames(datExpr)), datExpr)
rownames(datExpr_tmp)[1] <- "gene"
write.table(datExpr_tmp,file=paste0(datadir, "Expr_WholeBrain.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")

# read in MetaData
datMeta <- read.table(paste0(datadir, "Brain-AllTissues.samples_removed.regressed.covars.txt"), header=TRUE, row.names=1)
# fix sample names
sample_names <- gsub("-",".",rownames(datMeta))
rownames(datMeta) <- sample_names

# total brain without cerebellum network
samples <- rownames(datMeta[!(datMeta[,"tissue"] %in% c("CBH","CBL")),])
datExpr_subset <- datExpr[,samples]
# fix line 1 with sample names
datExpr_tmp <- rbind(c(colnames(datExpr_subset)), datExpr_subset)
rownames(datExpr_tmp)[1] <- "gene"
write.table(datExpr_tmp,file=paste0(datadir, "Expr_WholeBrainNoCB.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")

# Subset to each specific tissue
tissues <- unique(datMeta$tissue)
for (i in 1:length(tissues)) {
	tissue <- tissues[i]
	samples <- rownames(datMeta[datMeta[,"tissue"]==tissue,])
	datExpr_subset <- datExpr[,samples]
	# fix line 1 with sample names
	datExpr_tmp <- rbind(c(colnames(datExpr_subset)), datExpr_subset)
	rownames(datExpr_tmp)[1] <- "gene"
	write.table(datExpr_tmp,file=paste0(datadir, "Expr_",tissue,".txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
}