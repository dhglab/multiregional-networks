CRAN_REQ = c('zipfR', 'iterators', 'itertools', 'argparse')
BIOC_REQ = c('biomaRt', 'org.Hs.eg.db', 'topGO', 'clusterProfiler', 'DOSE', 'EnsDb.Hsapiens.v75', 
             'TxDb.Hsapiens.UCSC.hg19.knownGene', 'ChIPpeakAnno', 'WGCNA')
install.packages(CRAN_REQ, repos='http://cran.stat.ucla.edu')
source('http://bioconductor.org/biocLite.R')
biocLite(BIOC_REQ)
