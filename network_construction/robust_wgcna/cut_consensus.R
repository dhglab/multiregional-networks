library(WGCNA)
code.dir <- Sys.getenv('INTEGRATED_PIPELINE_DIR', '.')
source(paste(code.dir, 'integrated_pipeline_utils.R', sep='/'))

N_MODULES <- 6
N_BOOTSTRAPS <- 15
args <- commandArgs(trailingOnly=T)

consensus.tom.rda <- args[1]
bootstrap.dir <- args[2]
dat.expr.file <- args[3]
out.rdata <- args[4]

check.cuts <- rev(1 - exp(seq(log(0.005), log(0.2), length.out=25)))

check.names <- function(a, b) {
  if ( length(a) != length(b) ) {
    F
  } else {
    all(a == b)
  }
}

load(consensus.tom.rda)
colnames(consensus.TOM) <- rownames(consensus.TOM)

dat.expr <- read.table(dat.expr.file, header=T, row.names=1)
if ( ! check.names(rownames(dat.expr), rownames(consensus.TOM)) ) {
  msg <- sprintf('Names do not match:')
  de.names <- paste(head(rownames(dat.expr)), collapse=', ')
  tm.names <- paste(head(rownames(consensus.TOM)), collapse=', ')
  msg <- sprintf('%s\nde: %s\ntom: %s\n\n', msg, de.names, tm.names)
  stop(msg)
}

tree = hclust(as.dist(1-consensus.TOM), method='average')
modules.by.cut <- lapply(check.cuts, function(cut.height) {
   ct.res <- cutreeHybrid(dendro=tree, distM=1-consensus.TOM, deepSplit=2, pamStage=F,
                           pamRespectsDendro=F, minClusterSize=50, cut.height)
   ct.res$colors <- labels2colors(ct.res$labels)
   module.genes <- list()
   for ( label in unique(ct.res$colors) ) {
     module.genes[[label]] <- rownames(consensus.TOM)[ct.res$colors == label]
     print(sprintf('%s size %d', label, sum(ct.res$colors == label)))
   }
   ct.res[['module.genes']] <- module.genes
   ct.res
})

duplicate.modules <- list()
prev.table <- NULL
i <- 1
for ( cr in modules.by.cut ) {
  tab <- sort(table(cr$colors))
  if ( is.null(prev.table) ) {
    prev.table <- tab
  } else {
    if ( length(tab) == length(prev.table) && all(tab == prev.table) ) {
      duplicate.modules[[1 + length(duplicate.modules)]] <- i
    } else {
      prev.table <- tab
    }
  }
}

print(sprintf('%d equivalent cuts', length(duplicate.modules)))
non.duplicates <- setdiff(1:length(check.cuts), unlist(duplicate.modules))

modules.by.cut <- modules.by.cut[non.duplicates]

MeanBootTom <- function(bootstrap.file) {
  bf <- bootstrap.file
  # objective function: TOM correlation of the top three modules
  print(bf)
  boot.tom <- wrap.symm(scan(bf))
  sapply(1:length(modules.by.cut), function(idx) {
    mod.genes <- modules.by.cut[[idx]][['module.genes']]
    if ( length(mod.genes) <= N_MODULES ) { # 3 + grey = 4
      0
    } else {
      module.sizes <- sapply(mod.genes, length)
      top.modules <- order(module.sizes, decreasing=T)[1:(1+N_MODULES)]
      tom.cors <- rep(0, N_MODULES)
      j <- 1
      for ( mod in top.modules ) {
        print(names(mod.genes)[mod])
        if ( names(mod.genes)[mod] == 'grey' ) { # grey
          next
        }
        mgenes <- mod.genes[[mod]]
        cons.mtom <- consensus.TOM[mgenes,mgenes]
        mgene.idx <- which(rownames(consensus.TOM) %in% mgenes)
        boot.mtom <- boot.tom[mgene.idx, mgene.idx]
        tom.cors[j] <- cor(cons.mtom[lower.tri(cons.mtom)], boot.mtom[lower.tri(boot.mtom)])
        j <- 1 + j
        if ( j > N_MODULES ) {
          break
        }
      }
      mean(tom.cors)
   }
 })
}

MinBootTom <- function(bootstrap.file) {
 bf <- bootstrap.file
 print(bf)
 boot.tom <- wrap.symm(scan(bf))
 sapply(1:length(modules.by.cut), function(idx) {
    mod.genes <- modules.by.cut[[idx]][['module.genes']]
    tom.cors <- rep(0, length(mod.genes) - 1)
    j <- 1
    for ( mod in names(mod.genes) ) {
      if ( mod == 'grey' ) {
        next
      }
      mgenes <- mod.genes[[mod]]
        cons.mtom <- consensus.TOM[mgenes,mgenes]
        mgene.idx <- which(rownames(consensus.TOM) %in% mgenes)
        boot.mtom <- boot.tom[mgene.idx, mgene.idx]
        tom.cors[j] <- cor(cons.mtom[lower.tri(cons.mtom)], boot.mtom[lower.tri(boot.mtom)])
        j <- 1 + j
    }
    matrix(min(tom.cors), nrow=1,ncol=1)
  })
}

objective_function <- MeanBootTom

boot.files <- sprintf('%s/%s', bootstrap.dir, dir(bootstrap.dir))
boot.objectives <- lapply(boot.files[1:N_BOOTSTRAPS], objective_function)

# boot.objectives is a list of mean cor-TOMs by cut
boot.objectives <- do.call(rbind, boot.objectives)  # now looks like nboot x ncut
print(boot.objectives)
mean.objectives <- apply(boot.objectives, 2, mean)  # 1 x ncut
print(mean.objectives)
best.cut <- which.max(mean.objectives)

print(sprintf('best.cut=%f (%f)', check.cuts[best.cut], mean.objectives[best.cut]))

modules <- modules.by.cut[[best.cut]]
modules[['cut.height']] <- check.cuts[best.cut]
modules[['module.colors']] <- modules[['colors']]
names(modules[['module.colors']]) <- rownames(consensus.TOM)

mcm.res <- WGCNA::mergeCloseModules(t(dat.expr), modules$module.colors, relabel=T, corFnc='bicor')
names(mcm.res$colors) <- rownames(consensus.TOM)
modules$original.module.colors <- modules$module.colors
colors.to.labels <- sprintf('M%d', 1:60)
names(colors.to.labels) <- WGCNA::labels2colors(1:60)
modules$module.colors <- mcm.res$colors
modules$module.labels <- colors.to.labels[mcm.res$colors]

save(modules, file=out.rdata) 
