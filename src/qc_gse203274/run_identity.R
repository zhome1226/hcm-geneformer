suppressPackageStartupMessages({library(SingleR);library(SingleCellExperiment);library(Matrix);library(BiocParallel)})
root <- path.expand('~/xiaoe-internship-data')
work <- file.path(root,'_processing/all_data_cleaning_20261111/GSE203274')
refdir <- file.path(root,'processed/SCP1303_review_20261111/identity')
normalize <- function(x) {
    x <- as(x,'dgCMatrix');dn <- dimnames(x)
    x <- x %*% Diagonal(x=10000/pmax(colSums(x),1));x@x <- log1p(x@x);dimnames(x) <- dn;x
}
read_reference <- function(prefix) {
    con <- gzfile(file.path(refdir,paste0(prefix,'_counts.mtx.gz')),'rt');on.exit(close(con))
    as(readMM(con),'dgCMatrix')
}
prediction_frame <- function(pred) {
    scores <- as.matrix(pred$scores);assigned <- match(pred$labels,colnames(scores))
    data.frame(cell_id=rownames(pred),reference_label=pred$labels,reference_pruned=pred$pruned.labels,
               reference_delta=scores[cbind(seq_len(nrow(scores)),assigned)]-apply(scores,1,median),
               reference_delta_next=pred$delta.next,
               reference_best_score=scores[cbind(seq_len(nrow(scores)),assigned)])
}
set.seed(20261111)
common <- readLines(file.path(work,'reference_common_genes.txt'))
genes <- read.delim(file.path(refdir,'reference_genes.tsv'))$gene_id
rmeta <- read.delim(file.path(refdir,'reference_metadata.tsv'))
cmeta <- read.delim(file.path(refdir,'calibration_metadata.tsv'))
stopifnot(length(common)>5000,all(common %in% genes),!anyDuplicated(common))
ref <- read_reference('reference');cal <- read_reference('calibration')
rownames(ref) <- rownames(cal) <- genes
colnames(ref) <- rmeta$reference_barcode;colnames(cal) <- cmeta$reference_barcode
trained <- trainSingleR(normalize(ref[common,,drop=FALSE]),labels=rmeta$broad_type,
                        de.method='wilcox',de.n=50,aggr.ref=TRUE,num.threads=6,BPPARAM=MulticoreParam(3))
prediction <- classifySingleR(normalize(cal[common,,drop=FALSE]),trained,num.threads=6)
calibration <- prediction_frame(prediction)
calibration$truth <- cmeta$broad_type[match(calibration$cell_id,cmeta$reference_barcode)]
write.table(calibration,file.path(work,'reference_calibration.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
saveRDS(trained,file.path(work,'reference_trained.rds'))
query <- as(readMM(file.path(work,'strict_reference_query.mtx')),'dgCMatrix')
cells <- readLines(file.path(work,'strict_query_cells.txt'))
stopifnot(nrow(query)==length(common),ncol(query)==length(cells))
rownames(query) <- common;colnames(query) <- cells
if(length(cells)) {
    prediction <- classifySingleR(normalize(query),trained,num.threads=6)
    result <- prediction_frame(prediction)
} else {
    result <- calibration[FALSE,setdiff(colnames(calibration),'truth')]
}
write.table(result,file.path(work,'identity.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
writeLines(paste('Identity classified:',length(cells),'strict cells; published pediatric nuclei; inherited author doublet removal; residual screening recorded separately'),file.path(work,'IDENTITY_SUCCESS.txt'))
