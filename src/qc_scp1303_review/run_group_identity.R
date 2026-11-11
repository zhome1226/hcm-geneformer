suppressPackageStartupMessages({library(SingleR); library(Matrix); library(BiocParallel)})
out <- commandArgs(trailingOnly=TRUE)[[1]]
refdir <- dirname(out)
norm <- function(x) {
  dn <- dimnames(x)
  x <- as(x, 'CsparseMatrix') %*% Diagonal(x=10000/pmax(Matrix::colSums(x),1))
  x@x <- log1p(x@x); dimnames(x)<-dn; x
}
read <- function(path) { con<-gzfile(path,'rt'); on.exit(close(con)); as(readMM(con),'CsparseMatrix') }
genes <- read.delim(file.path(refdir,'reference_genes.tsv'))
ref <- read(file.path(refdir,'reference_counts.mtx.gz'))
meta <- read.delim(file.path(refdir,'reference_metadata.tsv'))
rownames(ref)<-genes$gene_id; colnames(ref)<-meta$reference_barcode
trained <- trainSingleR(norm(ref), labels=meta$broad_type, de.method='wilcox', de.n=50,
  aggr.ref=TRUE, num.threads=8, BPPARAM=MulticoreParam(4))
for (prefix in c('calibration','query')) {
  x<-read(file.path(out,paste0(prefix,'_group_counts.mtx.gz')))
  m<-read.delim(file.path(out,paste0(prefix,'_group_metadata.tsv')))
  ids<-if(prefix=='query') read.delim(file.path(out,'query_group_genes.tsv'))$gene_id else genes$gene_id
  rownames(x)<-ids; colnames(x)<-m$profile_id
  x<-x[genes$gene_id,,drop=FALSE]
  pred<-classifySingleR(norm(x), trained=trained, num.threads=8)
  stopifnot(identical(rownames(pred),m$profile_id))
  result<-data.frame(profile_id=rownames(pred), label=pred$labels,
    pruned_label=pred$pruned.labels, delta_next=pred$delta.next,
    as.matrix(pred$scores), check.names=FALSE)
  write.table(result,file.path(out,paste0(prefix,'_group_predictions.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='')
}
writeLines('Group identity review completed',file.path(out,'GROUP_SUCCESS.txt'))
