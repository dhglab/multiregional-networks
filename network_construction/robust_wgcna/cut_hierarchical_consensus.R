library(WGCNA)
options(stringsAsFactors=F)
code.dir <- Sys.getenv('INTEGRATED_PIPELINE_DIR', '.')
source(paste(code.dir, 'integrated_pipeline_utils.R', sep='/'))

args <- commandArgs(trailingOnly=T)

parent.tom.rda <- args[1]
left.child.tom.rda <- args[2]
right.child.tom.rda <- args[3]
out.rdata <- args[4]

if ( length(args) > 4 ) {
  expression.files <- unlist(strsplit(args[5], ','))
  print('Will produce consensus merges')
} else {
  expression.files <- NULL
}

check.cuts <- rev(1 - exp(seq(log(0.005), log(0.1), length.out=25)))

parent.env <- new.env()
left.env <- new.env()
right.env <- new.env()

load(parent.tom.rda, env=parent.env)
load(left.child.tom.rda, env=left.env)
load(right.child.tom.rda, env=right.env)

colnames(parent.env[[names(parent.env)[1]]]) <- rownames(parent.env[[names(parent.env)[1]]])
colnames(left.env[[names(left.env)[1]]]) <- rownames(left.env[[names(left.env)[1]]])
colnames(right.env[[names(right.env)[1]]]) <- rownames(right.env[[names(right.env)[1]]])


ptom <- function() {
  parent.env[[names(parent.env)[1]]]
}

ltom <- function() {
  left.env[[names(left.env)[1]]]
}

rtom <- function() {
  right.env[[names(right.env)[1]]]
}

if ( ! all(rownames(ptom()) == rownames(ltom())) ) {
  stop('Parent != left')
}

if ( ! all(rownames(ptom()) == rownames(rtom())) ) {
  stop('Parent != right')
}

if ( ! all(rownames(ltom()) == rownames(rtom())) ) {
  stop('right != left')
}

tree = hclust(as.dist(1 - ptom()), method='average')

modules.by.cut <- lapply(check.cuts, function(cut.height) {
   ct.res <- cutreeHybrid(dendro=tree, distM=1-ptom(), deepSplit=2, pamStage=F,
                           pamRespectsDendro=F, minClusterSize=50, cut.height)
   ct.res$colors <- labels2colors(ct.res$labels)
   module.genes <- list()
   for ( label in unique(ct.res$colors) ) {
     module.genes[[label]] <- rownames(ptom())[ct.res$colors == label]
     #print(sprintf('%s size %d', label, sum(ct.res$colors == label)))
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

N_MODULES=10
MeanTom <- function() {
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
        if ( names(mod.genes)[mod] == 'grey' ) { # grey
          next
        }
        mgenes <- mod.genes[[mod]]
        left.tom <- ltom()[mgenes, mgenes]
        right.tom <- rtom()[mgenes, mgenes]
        tom.cors[j] <- cor(left.tom[lower.tri(left.tom)], right.tom[lower.tri(right.tom)])
        j <- 1 + j
        if ( j > N_MODULES ) {
          break
        }
      }
      mean(tom.cors)
   }
 })
}

MinTom <- function() {
 sapply(1:length(modules.by.cut), function(idx) {
    mod.genes <- modules.by.cut[[idx]][['module.genes']]
    tom.cors <- rep(0, length(mod.genes) - 1)
    j <- 1
    for ( mod in names(mod.genes) ) {
      if ( mod == 'grey' ) {
        next
      }
        mgenes <- mod.genes[[mod]]
        left.mtom <- ltom()[mgenes, mgenes]
        right.mtom <- rtom()[mgenes, mgenes]
        tom.cors[j] <- cor(left.mtom[lower.tri(left.mtom)], right.mtom[lower.tri(right.mtom)])
        j <- 1 + j
    }
    weights <- sapply(tom.cors, function(tc) {
      if ( tc < 0 ) {
        3
      } else {
        1
      }
    })
    value <- sum(weights * tom.cors)/sum(weights)
    matrix(value, nrow=1,ncol=1)
  })
}

objective_function <- MeanTom

objective <- objective_function()

# boot.objectives is a list of mean cor-TOMs by cut
print(objective)
best.cut <- which.max(objective)

print(sprintf('best.cut=%f', check.cuts[best.cut]))

modules <- modules.by.cut[[best.cut]]
modules[['cut.height']] <- check.cuts[best.cut]
modules[['module.colors']] <- modules[['colors']]
names(modules[['module.colors']]) <- rownames(ptom())

if ( any(is.na(modules$module.colors)) ) {
  stop('Initial NA values in module colors?')
}

if ( ! is.null(expression.files) ) {
  merged.modules <- lapply(expression.files, function(expr.file) {
    print(sprintf("Merging %s", expr.file))
    expr.data <- read.table(expr.file, header=T, row.names=1)
    mods.merged <- WGCNA::mergeCloseModules(t(expr.data), modules$module.colors, relabel=T, corFnc='bicor')
    orig.colors <- unique(modules$module.colors)
    color.map <- sapply(orig.colors, function(x) { mods.merged$colors[modules$module.colors == x][1] })
    names(color.map) <- orig.colors
    lapply(color.map, function(dest) { names(color.map)[color.map == dest] })
  })
  mod.colors <- unique(modules$module.colors)
  # count module co-clusters
  co.clusters <- matrix(0, nrow=length(mod.colors), ncol=length(mod.colors))
  rownames(co.clusters) <- colnames(co.clusters) <- mod.colors
  for ( idx1 in 1:(length(mod.colors)-1) ) {
     mod.1 <- mod.colors[idx1]
    for ( idx2 in (1+idx1):length(mod.colors) ) {
      mod.2 <- mod.colors[idx2]
      for ( merge in merged.modules ) {
        for ( mods.together in merge ) {
          if ( mod.1 %in% mods.together && mod.2 %in% mods.together ) {
            co.clusters[idx1, idx2] <- 1 + co.clusters[idx1, idx2]
            co.clusters[idx2, idx1] <- co.clusters[idx1, idx2]
            break
          }
        }
      }
    }
  }
  co.clusters <- co.clusters/length(merged.modules)
  co.clusters['grey',] <- 0
  co.clusters[,'grey'] <- 0
  diag(co.clusters) <- 1
  options(digits=2)
  print(co.clusters)
  # co.clusters is now the co-clustering similarity matrix, which we can cluster on
  co.clust.tree <- hclust(as.dist(1 - co.clusters), 'complete')
  merged.names <- mod.colors
  branches <- as.factor(WGCNA::moduleNumber(dendro=co.clust.tree, cutHeight=0.8, minSize=1))
  uniqueBranches <- levels(branches)
  nBranches <- nlevels(branches)
  numOnBranch <- table(branches)

  for ( branch in 1:nBranches ) {
    if ( numOnBranch[branch] > 1 ) {
      modsOnThisBranch <- names(branches)[branches == uniqueBranches[branch]]
      for ( j in 2:length(modsOnThisBranch) ) {
        any.merged <- T
        merged.names[merged.names == modsOnThisBranch[j]] <- modsOnThisBranch[1]
      }
    }
  }

  merged.map <- merged.names
  names(merged.map) <- mod.colors
  print('merged.map:')
  print(merged.map)
  # now re-label by size
  orig.genes.in.module <- sapply(mod.colors, function(nm) { sum(modules$module.colors == nm) })
  names(orig.genes.in.module) <- mod.colors
  print('orig.genes')
  print(orig.genes.in.module)
  genes.in.module <- sapply(merged.names, function(nm) { sum(orig.genes.in.module[names(merged.map)[merged.map == nm]]) })
  names(genes.in.module) <- merged.names
  print('new.genes')
  print(genes.in.module)
  gim.uq <- unique(genes.in.module)
  gim.ranks <- rank(-gim.uq)
  module.ranks <- sapply(genes.in.module, function(x) { gim.ranks[gim.uq == x][1] })
  names(module.ranks) <- merged.names
  print('module.ranks')
  print(module.ranks)
  grey.num <- module.ranks['grey']
  module.ranks['grey'] <- 0
  module.ranks[module.ranks > grey.num] <- module.ranks[module.ranks > grey.num] - 1
  print(orig.genes.in.module)
  print(module.ranks)
  new.names <- sprintf('M%d', module.ranks)
  new.colors <- WGCNA::labels2colors(module.ranks)
  merged.names <- new.colors
  merged.map <- module.ranks
  names(merged.map) <- mod.colors
  modules$orig.labels <- modules$labels
  modules$orig.colors <- modules$colors
  modules$labels <- module.ranks[modules$colors]
  names(modules$labels) <- names(modules$colors)
  merged.map <- new.colors
  names(merged.map) <- mod.colors
  print('merged.map updated')
  print(merged.map)
  if ( any(is.na(modules$orig.colors)) ) {
    stop('NAs in original colors?')
  }
  if ( any(! modules$orig.colors %in% names(merged.map)) ) {
    print('Missing colors:')
    print(unique(modules$orig.colors[! modules$orig.colors %in% names(merged.map)]))
  } else {
    print('merged.map contains all original colors')
  }
  modules$colors <- merged.map[modules$orig.colors]
  names(modules$colors) <- names(modules$orig.colors)
  mod2lab <- 0:50
  names(mod2lab) <- WGCNA::labels2colors(mod2lab)
  modules$labels <- mod2lab[modules$colors]
  names(modules$labels) <- names(modules$colors)
}

if ( any(is.na(modules$colors)) || any(is.na(modules$labels)) ) {
  stop('NA values in final modules..?')
}
 

save(modules, file=out.rdata) 
