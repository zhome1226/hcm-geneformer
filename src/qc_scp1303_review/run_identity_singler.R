suppressPackageStartupMessages({library(SingleR); library(Matrix); library(BiocParallel)})
args <- commandArgs(trailingOnly=TRUE)
out <- args[[1]]
querydir <- dirname(out)
set.seed(20261111)
normalize <- function(x) {
  saved_names <- dimnames(x)
  x <- as(x, 'dgCMatrix')
  x <- x %*% Diagonal(x=10000/pmax(Matrix::colSums(x), 1))
  x@x <- log1p(x@x)
  dimnames(x) <- saved_names
  x
}
read_counts <- function(name) {
  con <- gzfile(file.path(out, paste0(name, '_counts.mtx.gz')), 'rt')
  on.exit(close(con))
  as(readMM(con), 'dgCMatrix')
}
genes <- read.delim(file.path(out, 'reference_genes.tsv'))
reference <- read_counts('reference')
calibration <- read_counts('calibration')
rmeta <- read.delim(file.path(out, 'reference_metadata.tsv'))
cmeta <- read.delim(file.path(out, 'calibration_metadata.tsv'))
stopifnot(nrow(reference)==nrow(genes), ncol(reference)==nrow(rmeta), ncol(calibration)==nrow(cmeta))
rownames(reference) <- rownames(calibration) <- genes$gene_id
colnames(reference) <- rmeta$reference_barcode
colnames(calibration) <- cmeta$reference_barcode
cat('Training SingleR on ', ncol(reference), ' reference nuclei\n', sep='')
trained <- trainSingleR(ref=normalize(reference), labels=rmeta$broad_type,
  de.method='wilcox', de.n=50, aggr.ref=TRUE,
  num.threads=8, BPPARAM=MulticoreParam(workers=4))
saveRDS(trained, file.path(out, 'singler_trained.rds'))
write_predictions <- function(pred, metadata, prefix) {
  scores <- as.matrix(pred$scores)
  assigned <- match(pred$labels, colnames(scores))
  delta <- scores[cbind(seq_len(nrow(scores)), assigned)] - apply(scores, 1, median)
  expected_barcodes <- if ('reference_barcode' %in% names(metadata)) metadata$reference_barcode else metadata$barcode
  stopifnot(identical(rownames(pred), expected_barcodes))
  df <- data.frame(barcode=rownames(pred), labels=pred$labels,
    pruned_labels=pred$pruned.labels, delta=delta, delta_next=pred$delta.next,
    best_score=scores[cbind(seq_len(nrow(scores)), assigned)])
  write.table(df, file.path(out, paste0(prefix, '_predictions.tsv')), sep='\t', row.names=FALSE, quote=FALSE, na='')
  write.table(data.frame(barcode=rownames(pred), scores, check.names=FALSE),
    file.path(out, paste0(prefix, '_scores.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
  saveRDS(pred, file.path(out, paste0(prefix, '_SingleR.rds')))
  df
}
cal <- classifySingleR(test=normalize(calibration), trained=trained, num.threads=8)
caldf <- write_predictions(cal, cmeta, 'calibration')
caldf$truth <- cmeta$broad_type[match(caldf$barcode, cmeta$reference_barcode)]
write.table(caldf, file.path(out, 'calibration_predictions.tsv'), sep='\t', row.names=FALSE, quote=FALSE, na='')
cat('Reference-donor holdout accuracy: ', mean(caldf$labels==caldf$truth), '\n', sep='')
rm(reference, calibration); gc()
cat('Reading query count-preserving matrix\n')
con <- gzfile(file.path(querydir, 'query_counts.mtx.gz'), 'rt')
query <- as(readMM(con), 'dgCMatrix'); close(con)
qgenes <- read.delim(file.path(querydir, 'query_genes.tsv'))
qmeta <- read.delim(file.path(querydir, 'query_metadata.tsv'), check.names=FALSE)
stopifnot(nrow(query)==nrow(qgenes), ncol(query)==nrow(qmeta), !anyDuplicated(qmeta$barcode))
rownames(query) <- qgenes$gene_id
colnames(query) <- qmeta$barcode
query <- query[genes$gene_id, , drop=FALSE]
pred <- classifySingleR(test=normalize(query), trained=trained, num.threads=8)
write_predictions(pred, qmeta, 'query')
capture.output(sessionInfo(), file=file.path(out, 'SINGLER_sessionInfo.txt'))
writeLines('SingleR reference calibration and query annotation completed', file.path(out, 'SINGLER_SUCCESS.txt'))
cat('SINGLER_IDENTITY_COMPLETE\n')
