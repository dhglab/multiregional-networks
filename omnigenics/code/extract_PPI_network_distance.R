### Combine and extract TFB data, turn these into networks, and compute distances upon them
options(stringsAsFactors=F)
args <- commandArgs(trailingOnly=T)

if ( length(args) == 0 || any(grepl('help|-h|--h|--help', args))) {
	stop('Usage: extract_PPI_network_distance.R [ppi_matrix] [output] [optional CORE1=gene1,gene2;CORE2=gene3,gene4')
}

library(biomaRt)
library(igraph)
mart <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

convert_to_ensembl <- function(gene.vec) {
	atts <- c('ensembl_gene_id', 'external_gene_name', 'external_gene_source', 'chromosome_name')
	gene.info <- getBM(atts, filters='external_gene_name', values=gene.vec, mart=mart)
	gene.info <- subset(gene.info, chromosome_name %in% c(1:22, 'X', 'Y', 'MT'))
	unique(gene.info$ensembl_gene_id)
}
convert_to_symbol <- function(gene.vec) {
  atts <- c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name')
  gene.info <- getBM(atts, filters='ensembl_gene_id', values=gene.vec, mart=mart)
  gene.info <- subset(gene.info, chromosome_name %in% c(1:22, 'X', 'Y', 'MT'))
  sym.vec <- gene.info$hgnc_symbol
  names(sym.vec) <- gene.info$ensembl_gene_id
  sym.vec[gene.vec]
}

library(igraph)
ppi_file <- args[1]
output_file <- args[2]
if ( length(args) >= 3 ) {
	raw_syndromic_genes <- args[3]
	sets = strsplit(raw_syndromic_genes, ';')[[1]]
	syndromic_gene_sets = list()
	for ( set in sets ) {
		setname <- strsplit(set, '=')[[1]][1]
		gene_str <- strsplit(set, '=')[[1]][2]
		genes <- strsplit(gene_str, ',')[[1]]
		syndromic_gene_sets[[setname]] <- convert_to_ensembl(genes)
	}
}
scale01 <- function(v) {
        (v - min(v))/(max(v)-min(v))
}

get_graph_distances <- function(graph, vertex.set, set.name) {
        # Extract mean and minimum distance from each vertex in the graph
        # to the given vertex set
        raw.distances <- igraph::distances(graph, to=intersect(vertex.set, V(graph)$name))
        raw.distances[is.infinite(raw.distances)] <- 2 * max(raw.distances[is.finite(raw.distances)])
        print(summary(raw.distances))
        collapsed.distances <- matrix(0, nrow=nrow(raw.distances), ncol=2)
        rownames(collapsed.distances) <- rownames(raw.distances)
        print(head(collapsed.distances))
        colnames(collapsed.distances) <- c(sprintf('%s.mean', set.name), sprintf('%s.min', set.name))
        collapsed.distances[,1] <- scale01(apply(raw.distances, 1, median))
        collapsed.distances[,2] <- scale01(apply(raw.distances, 1, function(x) { quantile(x, 3e-4) }))
        max.dist <- max(collapsed.distances[is.finite(collapsed.distances)])
        collapsed.distances[is.infinite(collapsed.distances)] <- 1 + max.dist
        collapsed.distances
}


ppi_dat <- as.matrix(read.table(ppi_file, header=T, row.names=1, sep='\t'))
ppi <- graph_from_adjacency_matrix(ppi_dat, mode='undirected')

degree.cent <- strength(ppi)
print(head(degree.cent))
central.genes <- names(degree.cent)[order(degree.cent, decreasing=T)[1:20]]
print(central.genes)

DISTANCES <- list()
DISTANCES[['central.genes.ppi']] <- get_graph_distances(ppi, central.genes, 'central.genes')


if ( length(args) >= 3 ) {
	for ( setname in names(syndromic_gene_sets) ) {
		print(sprintf(' .. .. .. %s .. .. ..', setname))
		set <- syndromic_gene_sets[[setname]]
		name <- sprintf('%s.ppi', setname)
		DISTANCES[[name]] <- get_graph_distances(ppi, set, name)
	}
}

print(' .. Ensuring proper order and writing .. ')

gene.order <- sort(V(ppi)$name)
for ( nm in names(DISTANCES) ) {
	DISTANCES[[nm]] <- DISTANCES[[nm]][gene.order,]
}

DISTANCES <- as.data.frame(do.call(cbind, DISTANCES))
ocn <- colnames(DISTANCES)
DISTANCES$gene <- rownames(DISTANCES)
DISTANCES$symbol <- convert_to_symbol(DISTANCES$gene)
rownames(DISTANCES) <- NULL
DISTANCES <- DISTANCES[,c('gene', 'symbol', ocn)]
print(head(DISTANCES))

write.table(DISTANCES, file=output_file, sep='\t', quote=F)
