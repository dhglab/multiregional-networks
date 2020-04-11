library('WGCNA')
code.dir <- Sys.getenv('INTEGRATED_PIPELINE_DIR', '.')
source(paste(code.dir, 'integrated_pipeline_utils.R', sep='/'))
# need data, power, sample.exclude, out.tom

MAX_GENES <- 20000
MAX_BLOCK <- 8000

args <- commandArgs(trailingOnly = T)
if ( is.null(args) || length(args) != 5 ) {
  stop('USAGE: BOOT_TOM.R [datafile] [power] [adjacency.type] [sample.to.exclude] [out.file]')
}
user.data <- args[1]
power <- args[2]
adj.type <- args[3]
seed.value <- args[4]
out.tom <- args[5]

if ( is.na(as.numeric(power)) ) {
  power <- as.numeric(scan(power, '%d'))
} else {
  power <- as.numeric(power)
}


if ( ! is.numeric(power) ) {
  power <- 6
  print('bad power, setting to 6')
}

set.seed(as.integer(seed.value))

data.matrix <- read.table(user.data)
if ( ncol(data.matrix) <= MAX_GENES ) {
  adjacency_mat = adjacency(data.matrix[sample.int(nrow(data.matrix), replace=T),], power=power, type=adj.type, corFnc='bicor')
  TOM = TOMsimilarity(adjacency_mat)
} else {
  # for speed, have to do this blockwise
  temp.name <- sprintf('%s.tom_tmp', out.tom)
  bwm <- blockwiseModules(data.matrix[sample.int(nrow(data.matrix), replace=T),], maxBlockSize=MAX_BLOCK, power=power,
                          networkType='signed', TOMType='unsigned', saveTOMs=T, saveTOMFileBase=temp.name, nThreads=4)
  TOM <- matrix(0, nrow=ncol(data.matrix), ncol=ncol(data.matrix))
  for ( block.idx in 1:length(bwm$TOMFiles) ) {
    tom.env <- new.env()
    tom.file <- bwm$TOMFiles[block.idx]
    tom.idx <- bwm$blockGenes[[block.idx]]
    print(tom.file)
    load(tom.file, env=tom.env)
    TOM[tom.idx, tom.idx] <- as.matrix(tom.env[['TOM']])
    system(sprintf('rm %s', tom.file))
  }
  rownames(TOM) <- colnames(TOM) <- colnames(data.matrix)
}

write.table(unwrap.symm(TOM), file=out.tom, quote=F, row.names=F, col.names=F)
