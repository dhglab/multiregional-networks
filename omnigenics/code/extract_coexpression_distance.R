# take co-expression modules and extract a matrix of core-gene distances
############
# The distances in question are:
## 1) 1 - module kME to each module
## 2) Minimum + mean path distance (adjacency) to max(2.5%, 5) hubs of each module
## 3) Minimum + mean path distance (adjacency) to syndromic gene sets
## 4) Minimum + mean path distance (1% + 1-NN)     to max(2.5%, 5) hubs of each modul
## 5) Minimum + mean path distance (1% + 1-NN)     to syndromic gene sets
############
#
# This produces a large matrix of
#
#          dist1   dist2   dist3  ....
#  gene1
#  gene2
#
####################################

# Rscript extract_coexpression_distance.R /u/project/geschwind/chartl/projects/network_grant/20171027_V5.5_eval/BRNCTX/BRNCTX.modules.mapped.txt /u/project/geschwind/chartl/projects/network_grant/20171027_V5.5_eval/BRNCTX/BRNCTX.expression.no_ext.txt test_pfc.distances.txt SCZ=ATPA13,FXYD,NRXN1,RB1CC1;ASD=FMRP,ANK2,SYNGAP1,CHD8,SHANK2,SHANK3,SCN2A


options(stringsAsFactors=F)
help.strs <- c('help', '--h', '--help', '-h', '-help')
args <- commandArgs(trailingOnly=T)
print(args)

if ( length(args) == 0 || any(help.strs %in% args ) ) {
	stop('Usage: extract_coexpression_distance.R [module file] [expression file] [output file] [optional NAME=GENE1,GENE2;NAME2=GENE3]')
}



library(WGCNA)
library(biomaRt)
library(igraph)
library(cppRouting)
library(parallel)
mart <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

convert_to_ensembl <- function(gene.vec) {
	atts <- c('ensembl_gene_id', 'external_gene_name', 'external_gene_source', 'chromosome_name')
	gene.info <- getBM(atts, filters='external_gene_name', values=gene.vec, mart=mart)
	gene.info <- subset(gene.info, chromosome_name %in% c(1:22, 'X', 'Y', 'MT'))
	unique(gene.info$ensembl_gene_id)
}

convert_to_hgnc <- function(gene.vec) {
	atts <- c('ensembl_gene_id', 'hgnc_symbol')
	gene.info <- getBM(atts, filters='ensembl_gene_id', values=gene.vec, mart=mart)
	e2h <- gene.info$hgnc_symbol
	names(e2h) <- gene.info$ensembl_gene_id
	e2h[gene.vec]
}

get_graph_distances <- function(graph, vertex.set, set.name) {
	# Extract mean and minimum distance from each vertex in the graph
	# to the given vertex set
	raw.distances <- igraph::distances(graph, to=intersect(vertex.set, V(graph)$name))
	collapsed.distances <- matrix(0, nrow=nrow(raw.distances), ncol=2)
	rownames(collapsed.distances) <- rownames(raw.distances)
	colnames(collapsed.distances) <- c(sprintf('%s.mean', set.name), sprintf('%s.min', set.name))
	collapsed.distances[,1] <- apply(raw.distances, 1, mean)
	collapsed.distances[,2] <- apply(raw.distances, 1, min)
	max.dist <- max(collapsed.distances[is.finite(collapsed.distances)])
	collapsed.distances[is.infinite(collapsed.distances)] <- 1 + max.dist
	collapsed.distances
}

get_matrix_distances <- function(dmat, vertex.set, set.name) {
	# Extract mean and minimum distance from each vertex in the graph
	# to the given vertex set
	vertex.set <- intersect(vertex.set, colnames(dmat))
	collapsed.distances <- matrix(0, nrow=nrow(dmat), ncol=2)
	rownames(collapsed.distances) <- rownames(dmat)
	colnames(collapsed.distances) <- c(sprintf('%s.mean', set.name), sprintf('%s.min', set.name))
	collapsed.distances[,1] <- apply(dmat[,vertex.set], 1, mean)
	collapsed.distances[,2] <- apply(dmat[,vertex.set], 1, min)
	collapsed.distances
}

mod_file <- args[1]
expression_file <- args[2]
output_file <- args[3]
if ( file.exists(output_file) ) {
  stop('Output file exists')
}

if ( length(args) > 3 ) {
	raw_syndromic_genes <- args[4]
	sets = strsplit(raw_syndromic_genes, ';')[[1]]
	syndromic_gene_sets = list()
	for ( set in sets ) {
		setname <- strsplit(set, '=')[[1]][1]
		gene_str <- strsplit(set, '=')[[1]][2]
		genes <- strsplit(gene_str, ',')[[1]]
		if ( any(! grepl('^ENSG', genes)) ) {
			ensg.genes <- genes[grepl('^ENSG', genes)]
			conv.genes <- convert_to_ensembl(genes[! grepl('^ENSG', genes)])
			genes <- c(ensg.genes, conv.genes)
		}
		syndromic_gene_sets[[setname]] <- genes
	}
}

mod.df <- read.table(mod_file, sep='\t', header=T)
mod.vec <- mod.df$module
names(mod.vec) <- mod.df$gene
expression <- read.table(expression_file, sep='\t', header=T)
overlap.genes <- intersect(mod.df$gene, rownames(expression))
if ( length(overlap.genes) < 5000 ) {
	stop('Fewer than 5000 genes overlap between expression and modules')
}
expression <- expression[overlap.genes,]
mod.vec <- mod.vec[overlap.genes]

print(' .. Computing kME and kWithin ..')
eigengenes <- moduleEigengenes(t(expression), mod.vec)
kME <- signedKME(t(expression), eigengenes$eigengenes)
kWithin <- intramodularConnectivity.fromExpr(t(expression), mod.vec, 
	networkType='signed', power=10)$kWithin

print(head(kWithin))
names(kWithin) <- rownames(expression)

DISTANCES <- list()
DISTANCES[['kME']] <- 1 - kME

mods_nogrey <- unique(mod.df$module)
mods_nogrey <- mods_nogrey[mods_nogrey != 'grey']

hub.genes <- lapply(mods_nogrey, function(mod) {
	mod.genes <- names(mod.vec)[mod.vec == mod]
	n <- as.integer(max(0.025*length(mod.genes), 5))
	mod.kwithin <- kWithin[mod.genes]
	names(mod.kwithin)[order(mod.kwithin, decreasing=T)[1:n]]
})
names(hub.genes) <- mods_nogrey

print(' .. Building adjacency matrix .. ')


raw_dist <- acos(bicor(t(expression)))


print(' .. Computing spherical distances .. ')

j <- 1
for ( mod in mods_nogrey ) {
	print(sprintf(' .. .. .. %s (%d/%d) .. .. ..', mod, j, length(mods_nogrey)))
	hubs <- hub.genes[[mod]]
	name <- sprintf('%s.spherical', mod)
	DISTANCES[[name]] <- get_matrix_distances(raw_dist, hubs, name)
	j <- j + 1
}

if ( length(args) > 3 ) {
	for ( setname in names(syndromic_gene_sets) ) {
		print(sprintf(' .. .. .. %s .. .. ..', setname))
		set <- syndromic_gene_sets[[setname]]
		name <- sprintf('%s.spherical', setname)
		DISTANCES[[name]] <- get_matrix_distances(raw_dist, set, name)
	}
}

print(' .. Building eps=.25% + 1-NN graph .. ')
nv <- nrow(raw_dist)
ne <- nv*(nv-1)/2 + nv
p_cutoff <- (0.025 * ne + nv)/ne
cutoff <- quantile(raw_dist, p_cutoff) # the closest edges
adj_1e1nn <- 1 * (raw_dist <= cutoff)
to.add <- list()
diag(adj_1e1nn) <- 0
for ( i in 1:nrow(adj_1e1nn) ) {
	degree <- sum(adj_1e1nn[i,])
	if ( degree == 0 ) {
		# connect the node to its nearest neighbor
		distances <- raw_dist[i,]
		distances[i] <- 1
		j <- which.min(distances)
		to.add[[1 + length(to.add)]] <- c(i, j)
	}
}
for ( ij in to.add ) {
	adj_1e1nn[ij[1], ij[2]] = 1
	adj_1e1nn[ij[2], ij[1]] = 1
}
adj_1e1nn <- graph_from_adjacency_matrix(adj_1e1nn, mode='undirected', weighted='blah')

print(' .. Computing sparse-neighbor-based distances .. ')

j <- 1
for ( mod in mods_nogrey ) {
	print(sprintf(' .. .. .. %s (%d/%d) .. .. ..', mod, j, length(mods_nogrey)))
	hubs <- hub.genes[[mod]]
	name <- sprintf('%s.sparse_1nn', mod)
	DISTANCES[[name]] <- get_graph_distances(adj_1e1nn, hubs, name)
	j <- j + 1
}

if ( length(args) > 3 ) {
	for ( setname in names(syndromic_gene_sets) ) {
		print(sprintf(' .. .. .. %s .. .. ..', setname))
		set <- syndromic_gene_sets[[setname]]
		name <- sprintf('%s.sparse_1nn', setname)
		DISTANCES[[name]] <- get_graph_distances(adj_1e1nn, hubs, name)
	}
}

#save(list=ls(), file='debug.Rda')
print(' .. Ensuring proper order and writing .. ')

gene.order <- names(mod.vec)
for ( nm in names(DISTANCES) ) {
	DISTANCES[[nm]] <- DISTANCES[[nm]][gene.order,]
}

ocn <- do.call(c, lapply(DISTANCES, colnames))
print(ocn)
DISTANCES <- do.call(cbind, DISTANCES)
colnames(DISTANCES) <- as.character(ocn)
DISTANCES <- data.frame(DISTANCES, check.names=F)
DISTANCES$gene <- gene.order
DISTANCES$symbol <- convert_to_hgnc(DISTANCES$gene)
rownames(DISTANCES) <- NULL
n <- ncol(DISTANCES)
DISTANCES <- DISTANCES[,c('gene', 'symbol', colnames(DISTANCES)[-c(n-1,n)])]
colnames(DISTANCES) <- sub('kME', 'kME.', colnames(DISTANCES))


write.table(DISTANCES, file=output_file, sep='\t', quote=F)
save(raw_dist, adj_1e1nn, file=sprintf('%s.graphs.Rdata', output_file))



