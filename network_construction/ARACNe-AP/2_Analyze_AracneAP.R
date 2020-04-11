# Script to analyze networks created by ARACNE-AP
# Gokul Ramaswami 10-12-2016
# Modeled mainly off a script from Chris Hartl: 20160819_analyze_ARACNe.inprogress.R

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)
# Load libraries
source("../../BrainGenesetAnnotation/R/utils/GO/GOUtils.R") # Chris Hartl's cell type and GO enrichments
library('igraph')
library('WGCNA')
library('ggplot2')

########## FUNCTIONS ######################
# Function to wrap text for barplot axis labels
wrapText <- function(x, len) { 
	sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
###########################################

# get arguments from command lines
argv <- commandArgs(trailingOnly = TRUE)

# set up tissue
Tissue <- "WholeBrain" #CURRENT TISSUE
if(length(argv) > 0) {
  Tissue = argv[1]
}

# data paths
datadir <- paste0("../processed_data/",Tissue,"/")
inputdatadir <- "../data/"
outputdir <- paste0("../figures/",Tissue,"/")

# read in the network
MIN_MI = 0.5
meanAR <- read.table(paste0(datadir,'network.txt'), header=T)

# prune some edges
meanAR <- subset(meanAR, MI >= MIN_MI)
meanAR.RegulatorNetwork <- subset(meanAR, Target %in% Regulator)
# remove regulator --> regulator edges
meanAR <- subset(meanAR, ! Target %in% Regulator)

# this is a bipartite graph. Sort it by the names.

meanAR <- meanAR[order(meanAR[,'Target']),]
meanAR <- meanAR[order(meanAR[,'Regulator']),]

# note that this is a directed graph; but only TF -> Targets, so the TFs are communities; make it undirected for finding higher modularity.
meanBipartiteGraph <- graph_from_edgelist(as.matrix(meanAR[,c('Regulator', 'Target')]), directed=F) # create graph from list of edges

# Three different clustering methods
meanWalktrapModules <- cluster_walktrap(meanBipartiteGraph, steps=4) # random walk
meanGreedyModules <- cluster_fast_greedy(simplify(meanBipartiteGraph)) # optimize modularity score
meanEigenModules <- cluster_leading_eigen(meanBipartiteGraph) # leading nonnegative eigenvector

# read in expression data and metadata
datExpr <- read.table(paste0(inputdatadir, "Expr_",Tissue,".txt"), header=TRUE, row.names=1)
datMeta <- read.table(paste0(inputdatadir, "Brain-AllTissues.samples_removed.regressed.covars.txt"), header=TRUE, row.names=1)
# fix sample names in datMeta
sample_names <- gsub("-",".",rownames(datMeta))
rownames(datMeta) <- sample_names
# Subset to samples in datExpr
datMeta = datMeta[colnames(datExpr),]

if (!file.exists(paste0(datadir,Tissue,"_MEs.RData"))) {

	# walktrap eigengenes
	meanWalktrapModules.Membership <- as.character(membership(meanWalktrapModules))
	names(meanWalktrapModules.Membership) <- names(membership(meanWalktrapModules))
	# remove small modules
	bad.WT <- names(table(meanWalktrapModules.Membership))[table(meanWalktrapModules.Membership) < 50] # Min module size of 50
	meanWalktrapModules.Membership[meanWalktrapModules.Membership %in% bad.WT] <- "0"
	meanWalktrapModules.Membership <- factor(meanWalktrapModules.Membership)
	# get mean expression level of module constituents
	datMean.WT <- datExpr[rownames(datExpr) %in% names(meanWalktrapModules.Membership),]
	datMean.WT <- datMean.WT[order(rownames(datMean.WT)),]
	meanWalktrapModules.Membership <- meanWalktrapModules.Membership[order(names(meanWalktrapModules.Membership))]
	# calculate module eigengenes
	meanWalktrapModules.eigengenes <- moduleEigengenes(t(datMean.WT), meanWalktrapModules.Membership, excludeGrey=T)

	## eigenmodule eigengenes
	meanEigenModules.Membership <- as.character(membership(meanEigenModules))
	names(meanEigenModules.Membership) <- names(membership(meanEigenModules))
	bad.EM <- names(table(meanEigenModules.Membership))[table(meanEigenModules.Membership) < 50]
	meanEigenModules.Membership[meanEigenModules.Membership %in% bad.EM] <- "0"
	meanEigenModules.Membership <- factor(meanEigenModules.Membership)
	datMean.EM <- datExpr[rownames(datExpr) %in% names(meanEigenModules.Membership),]
	datMean.EM <- datMean.EM[order(rownames(datMean.EM)),]
	meanEigenModules.Membership <- meanEigenModules.Membership[order(names(meanEigenModules.Membership))]
	meanEigenModules.eigengenes <- moduleEigengenes(t(datMean.EM), meanEigenModules.Membership, excludeGrey=T)


	## greedy module eigengenes
	meanGreedyModules.Membership <- as.character(membership(meanGreedyModules))
	names(meanGreedyModules.Membership) <- names(membership(meanGreedyModules))
	bad.GR <- names(table(meanGreedyModules.Membership))[table(meanGreedyModules.Membership) < 50]
	meanGreedyModules.Membership[meanGreedyModules.Membership %in% bad.GR] <- "0"
	meanGreedyModules.Membership <- factor(meanGreedyModules.Membership)
	meanGreedyModules.Membership <- meanGreedyModules.Membership[order(names(meanGreedyModules.Membership))]
	meanGreedyModules.eigengenes <- moduleEigengenes(t(datMean.EM), meanGreedyModules.Membership, excludeGrey=T)
	
	save(meanWalktrapModules.Membership, meanWalktrapModules.eigengenes, meanEigenModules.Membership, meanEigenModules.eigengenes, meanGreedyModules.Membership, meanGreedyModules.eigengenes, file=paste0(datadir,Tissue,"_MEs.RData"))
} else {
	load(paste0(datadir,Tissue,"_MEs.RData"))
}

## ME0 is the GREY MODULE!!

## make some plots (ME vs tissue type)
TISSUE_STR <- sapply(colnames(datExpr), function(x) { strsplit(x, '\\.')[[1]][3] })

# Walktrap
pdf(paste0(outputdir,Tissue,"_Walktrap_Tissue.pdf"))
for ( wtEig.str in sort(colnames(meanWalktrapModules.eigengenes$eigengenes)) ) {
  wtEig <- meanWalktrapModules.eigengenes$eigengenes[,wtEig.str]
  df <- data.frame(eigenPC=wtEig, tissue=TISSUE_STR)
  n <- paste('(n=', table(meanWalktrapModules.Membership)[substr(wtEig.str, 3, nchar(wtEig.str))], ')', sep='')
  x <- ggplot(df) + geom_boxplot(aes_string(x='tissue', y='eigenPC', fill='tissue')) + ggtitle(paste('Walktrap', wtEig.str, n))
  print(x)
}
dev.off()

# MeanEigen
pdf(paste0(outputdir,Tissue,"_MeanEigen_Tissue.pdf"))
for ( meEig.str in sort(colnames(meanEigenModules.eigengenes$eigengenes)) ) {
  meEig <- meanEigenModules.eigengenes$eigengenes[,meEig.str]
  df <- data.frame(eigenPC=meEig, tissue=TISSUE_STR)
  n <- paste('(n=', table(meanEigenModules.Membership)[substr(meEig.str, 3, nchar(meEig.str))], ')', sep='')
  x <- ggplot(df) + geom_boxplot(aes_string(x='tissue', y='eigenPC', fill='tissue')) + ggtitle(paste('MeanEigen', meEig.str, n))
  print(x)
}
dev.off()

# Greedy
pdf(paste0(outputdir,Tissue,"_Greedy_Tissue.pdf"))
for ( grEig.str in sort(colnames(meanGreedyModules.eigengenes$eigengenes)) ) {
  grEig <- meanGreedyModules.eigengenes$eigengenes[,grEig.str]
  df <- data.frame(eigenPC=grEig, tissue=TISSUE_STR)
  n <- paste('(n=', table(meanGreedyModules.Membership)[substr(grEig.str, 3, nchar(grEig.str))], ')', sep='')
  x <- ggplot(df) + geom_boxplot(aes_string(x='tissue', y='eigenPC', fill='tissue')) + ggtitle(paste('GreedyEigen', grEig.str, n))
  print(x)
}
dev.off()

# Plot ME vs traits (Age, Ethnicity, Gender, Race, SMGEBTCH, SMNABTCH, RIN)
methodNames <- c("Walktrap","Meaneigen","Greedy")
MEs_list <- list(meanWalktrapModules.eigengenes$eigengenes,meanEigenModules.eigengenes$eigengenes,meanGreedyModules.eigengenes$eigengenes)
for (j in 1:length(MEs_list)) {
	MEs <- MEs_list[[j]]
	Pmat <- matrix(NA,nrow=ncol(MEs),ncol=3)
	colnames(Pmat) = c("AGE","GENDER","RIN")
	rownames(Pmat) = colnames(MEs)
	Bmat <- matrix(NA,nrow=ncol(MEs),ncol=3)
	for (i in 1:ncol(MEs)) {
		# LINEAR MODEL 
		age_lm = lm(MEs[,i]~datMeta[,"AGE"])
		gender_lm = lm(MEs[,i]~as.factor(datMeta[,"GENDER"]))
		rin_lm = lm(MEs[,i]~datMeta[,"SMRIN"])
		Pmat[i,] = c(summary(age_lm)$coefficients[2,4],summary(gender_lm)$coefficients[2,4],summary(rin_lm)$coefficients[2,4])
		Bmat[i,] = c(summary(age_lm)$coefficients[2,1],summary(gender_lm)$coefficients[2,1],summary(rin_lm)$coefficients[2,1])
	}
	pdf(paste0(outputdir,Tissue,"_",methodNames[j],"_MEtraits.pdf"))
	## Plot the module-level summaries
	BFmat <- Pmat*(nrow(Pmat)-1)
	BFmat[BFmat>1] <- 1
	FDRmat <- BFmat ## USING BONFERRONI - comment this out if you don't want to!

	textFDRmat <- signif(FDRmat,2)
	dispMat <- -log10(Pmat)*sign(Bmat)
	dispMat[FDRmat > 0.05] <- 0
	dispMat[dispMat > 10] <- 10
	dispMat[dispMat < -10] <- -10
	rownames(dispMat) <- colnames(MEs)

	maxval <- 10
	labeledHeatmap(Matrix=dispMat,
               textMatrix=textFDRmat,
               yLabels=colnames(MEs),
               yColorLabels=TRUE,
               ySymbols=colnames(MEs),
               xLabels=colnames(Pmat),
               colors=blueWhiteRed(100),
               cex.lab.x=1,
               cex.lab.y=0.7,
               cex.text=0.5,
               zlim=c(-maxval,maxval),
               main="ME - trait associations")
	dev.off()
}


# loop over 3 grouping methods for GO and cell type enrichments
colors_list <- list(meanWalktrapModules.Membership,meanEigenModules.Membership,meanGreedyModules.Membership)

# Cell type enrichments
# loop over 3 grouping methods
for (i in 1:length(colors_list)) {
	pdf(paste0(outputdir,Tissue,"_",methodNames[i],"_CellTypeEnrich.pdf"))
	moduleColors <- colors_list[[i]]
	# GO enrichments for modules
	modules <- unique(moduleColors)
	for (module in 1:length(modules)) { 
		mod.ind <- 1 * (moduleColors==modules[module])
		mod.scores <- mod.ind
		names(mod.ind) <- names(moduleColors)
		names(mod.scores) <- names(moduleColors)
		mod.cell.type <- ModuleBrainCellEnrich(mod.ind, mod.scores)
		par(mar = c(5,7,4,2) + 0.1)
		barplot(-log10(mod.cell.type[,"fdr.score"]),horiz=TRUE,names.arg = 	wrapText(mod.cell.type[,"term"],20),xlab="-log10(pvalue)",cex.names=0.5,col="lightblue",las=2,main=paste0("M.",modules[module]))
		abline(v=-log10(0.05),col="red")
	}
	dev.off()
}

# GO enrichments
# loop over 3 grouping methods
for (i in 1:length(colors_list)) {
	pdf(paste0(outputdir,Tissue,"_",methodNames[i],"_GOenrich.pdf"))
	moduleColors <- colors_list[[i]]
	# GO enrichments for modules
	modules <- unique(moduleColors)
	for (module in 1:length(modules)) { 
		mod.ind <- 1 * (moduleColors==modules[module])
		mod.scores <- mod.ind
		names(mod.ind) <- names(moduleColors)
		names(mod.scores) <- names(moduleColors)
		mod.GO <- ModuleGOEnrich(mod.ind, mod.scores,min.genes=50,verbose=T)
		if (nrow(mod.GO) > 10) { # restrict to top 10 GO terms
			mod.GO = mod.GO[1:10,]
		}
		mod.GO.graph <- PlotGOEnrich(go.result  = mod.GO, 
                          p.cut      = 0.05, 
                          p.use      = 'p.ind', 
                          bar.length = '-log10(p.ind)',
                          order.by   = 'p.ind',
                          order.dec  = TRUE,
                          max.char   = 30,
                          fill.color = 'or.ind',
                          flip.axes  = TRUE,
                          min.terms  = 5,
                          title_label = paste0("M.",modules[module])) # I modified the source plotting code to include a title
        print(mod.GO.graph)
	}
	dev.off()
}

# KEGG enrichments
# loop over 3 grouping methods
for (i in 1:length(colors_list)) {
	pdf(paste0(outputdir,Tissue,"_",methodNames[i],"_KEGGenrich.pdf"))
	moduleColors <- colors_list[[i]]
	# GO enrichments for modules
	modules <- unique(moduleColors)
	for (module in 1:length(modules)) { 
		mod.ind <- 1 * (moduleColors==modules[module])
		mod.scores <- mod.ind
		names(mod.ind) <- names(moduleColors)
		names(mod.scores) <- names(moduleColors)
		mod.KEGG <- ModuleGOEnrich(mod.ind, mod.scores,min.genes=50,ontology='KEGG',verbose=T)
		if (nrow(mod.KEGG) > 10) { # restrict to top 10 KEGG terms
			mod.KEGG = mod.KEGG[1:10,]
		}
		mod.KEGG.graph <- PlotGOEnrich(go.result  = mod.KEGG, 
                          p.cut      = 0.05, 
                          p.use      = 'p.ind', 
                          bar.length = '-log10(p.ind)',
                          order.by   = 'p.ind',
                          order.dec  = TRUE,
                          max.char   = 30,
                          fill.color = 'or.ind',
                          flip.axes  = TRUE,
                          min.terms  = 5,
                          title_label = paste0("M.",modules[module])) # I modified the source plotting code to include a title
        print(mod.KEGG.graph)
	}
	dev.off()
}