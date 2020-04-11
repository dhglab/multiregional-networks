library('WGCNA')

code.dir <- Sys.getenv('INTEGRATED_PIPELINE_DIR', '.')
source(paste(code.dir, 'integrated_pipeline_utils.R', sep='/'))
# need data, power, sample.exclude, out.tom

args <- commandArgs(trailingOnly = T)
if ( is.null(args) || length(args) != 4 ) {
  stop('USAGE: LOO_TOM.R [datafile] [power] [sample.to.exclude] [out.file]')
}
user.data <- args[1]
power <- as.integer(args[2])
adj.type <- args[3]
sample.exclude <- args[4]
out.tom <- args[4]

if ( ! is.na(as.integer(sample.exclude)) ) {
  sample.exclude <- as.integer(sample.exclude)
}

data.matrix <- read.table(user.data)
adjacency_mat = adjacency(data.matrix[-sample.exclude,], power=power, type=adj.type)
TOM = TOMsimilarity(adjacency_mat)
write.table(unwrap.symm(TOM), file=out.tom, quote=F, row.names=F, col.names=F)
