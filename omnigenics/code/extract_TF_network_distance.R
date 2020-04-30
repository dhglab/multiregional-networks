### Combine and extract TFB data, turn these into networks, and compute distances upon them
options(stringsAsFactors=F)
args <- commandArgs(trailingOnly=T)

# Rscript extract_TF_network_distance.R /u/project/geschwind/chartl/projects/network_grant/20190326_omnigenics/data/regulatory_networks/hippocampus_adult.txt.gz,/u/project/geschwind/chartl/projects/network_grant/20190326_omnigenics/data/regulatory_networks/amygdala_adult.txt.gz test_TF.distance.txt "ASD=FMRP,ANK2,SHANK2,SHANK3;SCZ=ATP1A3,FXYD,NRXN1,RB1CC1"
if ( length(args) == 0 || any(grepl('help|-h|--h|--help', args))) {
	stop('Usage: extract_TF_network_distance.R [edge_list1,edge_list2] [output] [optional CORE1=gene1,gene2;CORE2=gene3,gene4')
}

library(biomaRt)
library(igraph)
mart <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

convert_to_ensembl <- function(gene.vec) {
	atts <- c('ensembl_gene_id', 'external_gene_name', 'external_gene_source', 'chromosome_name')
	gene.info <- getBM(atts, filters='external_gene_name', values=gene.vec, mart=mart)
	gene.info <- subset(gene.info, chromosome_name %in% c(1:22, 'X', 'Y', 'MT'))
	sapply(gene.vec, function(g) {
		gene.info[gene.info$external_gene_name == g,]$ensembl_gene_id[1]
	})
}

scale01 <- function(v) {
	(v - min(v))/(max(v)-min(v))
}

get_graph_distances <- function(graph, vertex.set, set.name) {
	# Extract mean and minimum distance from each vertex in the graph
	# to the given vertex set
	raw.distances <- igraph::distances(graph, to=intersect(vertex.set, V(graph)$name))
	collapsed.distances <- matrix(0, nrow=nrow(raw.distances), ncol=2)
	rownames(collapsed.distances) <- rownames(raw.distances)
	print(head(collapsed.distances))
	colnames(collapsed.distances) <- c(sprintf('%s.mean', set.name), sprintf('%s.min', set.name))
	collapsed.distances[,1] <- scale01(apply(raw.distances, 1, mean))
	collapsed.distances[,2] <- scale01(apply(raw.distances, 1, min))
	max.dist <- max(collapsed.distances[is.finite(collapsed.distances)])
	collapsed.distances[is.infinite(collapsed.distances)] <- 1 + max.dist
	collapsed.distances
}

library(igraph)

if ( length(args) > 2 ) {
	print(args)
	raw_syndromic_genes <- args[3]
	print(raw_syndromic_genes)
	sets = strsplit(raw_syndromic_genes, ';')[[1]]
	print(sets)
	syndromic_gene_sets = list()
	for ( set in sets ) {
		setname <- strsplit(set, '=')[[1]][1]
		gene_str <- strsplit(set, '=')[[1]][2]
		genes <- strsplit(gene_str, ',')[[1]]
		syndromic_gene_sets[[setname]] <- genes
	}
}



edge_files <- strsplit(args[1], ',')[[1]]
output_file <- args[2]

write.table(rnorm(100), file=output_file, sep='\t', quote=F)  # make sure it exists

edge.lists <- lapply(edge_files, function(fn) {
	read.table(fn, header=F)
})

for ( i in 1:length(edge.lists) ) {
	edge.lists[[i]]$key <- sprintf('%s-%s', edge.lists[[i]]$V1, edge.lists[[i]]$V2)
	colnames(edge.lists[[i]]) <- c('TF', 'gene', sprintf('score_%d', i), 'key')
}

edge.df <- Reduce(function(df1, df2){ merge(df1, df2, by=c('key', 'TF', 'gene'), all=T) }, edge.lists)
print(head(edge.df))
score.cols <- grepl('score', colnames(edge.df))
edge.df$mean_score <- apply(edge.df[,score.cols,drop=F], 1, mean, na.rm=T)
graph <- graph_from_data_frame(edge.df[,c('TF', 'gene', 'mean_score')], directed=F)# marbach does this
E(graph)$weight <- E(graph)$mean_score
print('Graph loaded, computing RWK')
M <- as_adj(graph, attr = 'mean_score', sparse = T) 
s <- sqrt(strength(graph, mode='all'))
D <- diag(1/s)
print(' .. (1)')
K <- D %*% M %*% D + diag(rep(1, nrow(M)))
print(' .. (2)')
K <- K %*% K  # (I+W)^2
print(' .. (3)')
K <- K %*% K  # (I+W)^4
print('done!')
rownames(K) <- colnames(K) <- rownames(M)
maxK <- max(K)
minK <- min(K)
kernel.graph <- graph_from_adjacency_matrix(K, mode='undirected', weighted=T, diag = F)
rm(graph, K, M, D)

# a natural choice of core genes are those with the highest centrality under
# the random walk kernel

print('degree centrality')
degree.cent <- strength(kernel.graph)
central.genes <- names(degree.cent)[order(degree.cent, decreasing=T)[1:25]]


print('Similarity -> Dissimilarity weights for distances')
E(kernel.graph)$sim.weight <- E(kernel.graph)$weight
E(kernel.graph)$weight <- 1 - (E(kernel.graph)$sim.weight - minK)/(maxK-minK)
# this is too expensive, so keep core genes and 5000 additional
all.genes <- Reduce(union, syndromic_gene_sets)
all.genes <- union(all.genes, read.table('omnigenics/data2/keep_syndromic_genes.txt', header=F)$V1)
rnd.genes <- sample(V(kernel.graph)$name, 4500)
all.genes <- union(all.genes, rnd.genes)
all.genes <- intersect(all.genes, V(kernel.graph)$name)

of <- strsplit(output_file, '/')[[1]]
of <- of[length(of)]
save(list=ls(), file=sprintf('debug.%s.Rda', of))

print(sprintf('For speed subsetting to %d/%d', length(all.genes), length(V(kernel.graph))))
kernel.graph <- subgraph(kernel.graph, all.genes)


DISTANCES <- list()
DISTANCES[['central.genes']] <- get_graph_distances(kernel.graph, central.genes, 'central.genes')


if ( length(args) > 2 ) {
	for ( setname in names(syndromic_gene_sets) ) {
		print(sprintf(' .. .. .. %s .. .. ..', setname))
		set <- syndromic_gene_sets[[setname]]
		name <- sprintf('%s.sparse_1nn', setname)
		DISTANCES[[name]] <- get_graph_distances(kernel.graph, set, name)
	}
}

print(' .. Ensuring proper order and writing .. ')

gene.order <-sort(Reduce(intersect,lapply(DISTANCES, rownames)))
for ( nm in names(DISTANCES) ) {
	DISTANCES[[nm]] <- DISTANCES[[nm]][gene.order,]
}

DISTANCES <- as.data.frame(do.call(cbind, DISTANCES))
ocn <- colnames(DISTANCES)
DISTANCES$symbol <- rownames(DISTANCES)
DISTANCES$gene <- convert_to_ensembl(DISTANCES$symbol)
rownames(DISTANCES) <- NULL
DISTANCES <- DISTANCES[,c('gene', 'symbol', ocn)]


write.table(DISTANCES, file=output_file, sep='\t', quote=F)



