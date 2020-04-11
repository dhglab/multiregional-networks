library('argparse')
library('WGCNA')

parser <- ArgumentParser()
parser$add_argument("-e", "--expression", type="character", help="The expression data file. Rows are genes, columns are samples.")
parser$add_argument("-o", "--output", type="character", help="The output directory", default="./")
parser$add_argument("-t", "--transpose", action="store_true", default=FALSE, help="the expression file has rows as samples, so transpose it")
parser$add_argument("--noHeader", action="store_true", default=FALSE, help="The expression data has no header")
args <- parser$parse_args()

if ( is.null(args$expression) ) {
  stop('Must specify expression data')
}

datExpr <- read.table(args$expression, header=!args$noHeader)
if ( rownames(datExpr)[1] == "1" ) {
  rownames(datExpr) <- datExpr[,1]
  datExpr <- datExpr[,-1]
}
if ( args$transpose ) {
  datExpr <- t(datExpr)
}

POWS <- 4:20
sft <- pickSoftThreshold(t(datExpr), powerVector=POWS, verbose=5, RsquaredCut=0.8, networkType='signed')
indepCriterion <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
outfile = paste(args$output, "WGCNA_power.pdf", sep='/')
pdf(outfile)
plot(POWS, indepCriterion, xlab='WGCNA power', ylab='Power Criterion (R^2 * sign of coefficient)', bty='n', col='red', type='o')
dev.off()

if ( any(indepCriterion > 0.8) ) {
  param.idx <- which(indepCriterion > 0.8)[1]
  param <- sft$fitIndices[param.idx,1]
} else {
  indepCriterion <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,4]
  if ( any(indepCriterion > 0.9) ) {
    param.idx <- which(indepCriterion > 0.9)[1]
    param <- sft$fitIndices[param.idx, 1]
  } else {
    param <- 20
  }
}

out_param = paste(args$output, 'power.param', sep='/')
cat(param, file=out_param)

print(warnings())
