suppressPackageStartupMessages({library(scDblFinder); library(SingleCellExperiment); library(Matrix); library(BiocParallel)})
out <- commandArgs(trailingOnly=TRUE)[[1]]
captures <- read.delim(file.path(out, 'capture_inputs.tsv'))
all_results <- BiocParallel::bplapply(seq_len(nrow(captures)), function(i) {
  capture <- captures$capture_id[[i]]
  directory <- file.path(out, capture)
  cells <- read.delim(file.path(directory, 'cells.tsv'))
  con <- gzfile(file.path(directory, 'counts.mtx.gz'), 'rt')
  counts <- as(readMM(con), 'CsparseMatrix'); close(con)
  colnames(counts) <- cells$barcode
  stopifnot(ncol(counts)==nrow(cells))
  set.seed(20261111+i)
  sce <- SingleCellExperiment(assays=list(counts=counts))
  sce <- scDblFinder(sce, dbr=0.01, dbr.sd=0.015, BPPARAM=SerialParam(), verbose=FALSE)
  result <- data.frame(barcode=colnames(sce), capture_id=capture,
    scDblFinder_score=sce$scDblFinder.score, scDblFinder_class=as.character(sce$scDblFinder.class))
  write.table(result, file.path(directory, 'predictions.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
  cat('Completed ', capture, '\n', sep='')
  result
}, BPPARAM=MulticoreParam(workers=4, RNGseed=20261111, progressbar=FALSE))
predictions <- do.call(rbind, all_results)
stopifnot(!anyDuplicated(predictions$barcode))
con <- gzfile(file.path(out, 'secondary_doublet_predictions.tsv.gz'), 'wt')
write.table(predictions, con, sep='\t', row.names=FALSE, quote=FALSE); close(con)
capture.output(sessionInfo(), file=file.path(out, 'R_sessionInfo.txt'))
writeLines('Secondary scDblFinder review completed', file.path(out, 'DOUBLET_SUCCESS.txt'))
