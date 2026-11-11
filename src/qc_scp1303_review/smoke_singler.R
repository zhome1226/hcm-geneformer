args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Matrix))
set.seed(20261111)
reference <- matrix(rpois(100 * 20, 2), nrow = 100)
reference[1:30, 1:10] <- reference[1:30, 1:10] + 30
reference[31:60, 11:20] <- reference[31:60, 11:20] + 30
rownames(reference) <- paste0('gene', seq_len(100))
labels <- rep(c('A', 'B'), each = 10)
query <- cbind(rowMeans(reference[, 1:10]), rowMeans(reference[, 11:20]))
rownames(query) <- rownames(reference)
colnames(query) <- c('synthetic_A', 'synthetic_B')
prediction <- SingleR(test = log1p(query), ref = log1p(reference), labels = labels)
stopifnot(identical(as.character(prediction$labels), c('A', 'B')))
write.table(as.data.frame(prediction[, c('labels', 'pruned.labels')]),
            file.path(out, 'singler_smoke_predictions.tsv'), sep = '\t', quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'R_sessionInfo.txt'))
writeLines(c('status=passed', paste0('SingleR=', packageVersion('SingleR')),
             paste0('SingleCellExperiment=', packageVersion('SingleCellExperiment')),
             paste0('R=', getRversion())), file.path(out, 'R_SMOKE_TEST.txt'))
cat('SingleR synthetic annotation smoke test passed\n')
