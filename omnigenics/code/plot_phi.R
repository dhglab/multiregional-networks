setwd('~/repos/multiregional-networks/')
rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)

core.mapping <- c(
  'DuWu'='NPDenovo',
  'extTADA_ASD'='extTADA.A',
  'extTADA_SCZ'='extTADA.S',
  'iHart'='iHart'
)

color.map <- c('extTADA.A'='#27d9d9',
               'extTADA.S'='#00c5ff',
               'NPDenovo'='#71a2ff',
               'iHart'='#de69d8',
               'Module'='#ff2e6d')


mod.pat <- c('WHOLE_BRAIN.M', '.spherical.min', '.spherical.mean',
             '.sparse_1nn', 'kME.', '.min', '.max', '.mean', 'BRNCTX.M',
             '.sparse_1nn.mean')
mod.tar <- c('BW-M', '', '', '', '', '', '', '', 'PFC-M', '')


map.core <- function(v) {
  x <- v
  for ( n in names(core.mapping) ) {
    x[grepl(n, x)] <- core.mapping[n]
  }
  x
}

map.dist <- function(v) {
  x <- map.core(v)
  mod.idx <- which(x == v)
  for ( i in 1:length(mod.pat) ) {
    x[mod.idx] <- gsub(mod.pat[i], mod.tar[i], x[mod.idx], fixed=T)
  }
  x
  
}

phi.dat <- c()
for ( phi.file in dir('omnigenics/outputs') ) {
  fn <- sprintf('omnigenics/outputs/%s', phi.file)
  if ( grepl('\\.phi\\.txt', fn) ) {
    nm <- strsplit(phi.file, '\\.')[[1]][1]
    dat <- read.table(fn, header=T, row.names = NULL)
    dat$network <- nm
    phi.dat <- rbind(phi.dat, dat)
  }
}

phi.dat <- phi.dat[order(phi.dat$phi, decreasing=T),]

phi.dat$core.name <- map.core(phi.dat$core)
phi.dat$dist.name <- map.dist(phi.dat$distance)
phi.dat$dist.type <- sapply(phi.dat$dist.name, function(n) {
  if ( grepl('M', n) ) {
    'Module'
  } else {
    n
  }
})
phi.dat$dist.type <- factor(phi.dat$dist.type, 
                                levels=c('extTADA.A', 'extTADA.S', 'NPDenovo', 'iHart', 'Module', 'grey'))

phi.dat$p.value <- sapply(1:nrow(phi.dat), function(i) {
  nG <- 10 * phi.dat$n_D[i]
  mat <- matrix(c(
    nG - phi.dat$n_D[i], phi.dat$n_D[i],
    phi.dat$n_core[i] - phi.dat$overlap[i], phi.dat$overlap[i]
  ), nrow=2
  )
  fisher.test(mat, alternative='greater')$p.value
})

valid.phi <- subset(phi.dat, n_core > 14 & grepl('protein', core))
valid.phi$network <- factor(valid.phi$network, levels=c('WHOLE_BRAIN', 'InWeb', 'RC_Neuron'))


pdf('phi_plots.pdf', height=4, width=8)
for ( cn in unique(valid.phi$core.name) ) {
  ddf <- subset(valid.phi,  core.name == cn & dist.type %in% names(color.map))
  # add 1 dummy record to ensure plots
  for ( nn in unique(valid.phi$network) ) {
    rec <- ddf[nrow(ddf),]
    rec$network <- nn
    rec$phi <- 0
    ddf <- rbind(ddf, rec)
  }
  gp <- ggplot(ddf)
  gp <- gp + geom_bar(aes(x=network, fill=dist.type, y=phi), stat='identity', position=position_dodge2(preserve='single'))
  gp <- gp + theme_bw() + ggtitle(cn)
  gp <- gp + scale_fill_manual(values = color.map)
  print(gp)
}
dev.off()

