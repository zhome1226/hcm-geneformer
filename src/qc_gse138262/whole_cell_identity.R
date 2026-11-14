suppressPackageStartupMessages({library(SingleR);library(Matrix);library(BiocParallel)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
work <- args[[1]]
refdir <- '/home/zhome/xiaoe-internship-data/processed/SCP1303_review_20261111/identity'
normalize <- function(x) {
  x <- as(x,'dgCMatrix'); dn <- dimnames(x)
  x <- x %*% Diagonal(x=10000/pmax(colSums(x),1)); x@x <- log1p(x@x); dimnames(x) <- dn; x
}
read_reference <- function(prefix) {
  con <- gzfile(file.path(refdir,paste0(prefix,'_counts.mtx.gz')),'rt'); on.exit(close(con))
  as(readMM(con),'dgCMatrix')
}
prediction_frame <- function(pred) {
  scores <- as.matrix(pred$scores); assigned <- match(pred$labels,colnames(scores))
  data.frame(cell_id=rownames(pred),reference_label=pred$labels,reference_pruned=pred$pruned.labels,
    reference_delta=scores[cbind(seq_len(nrow(scores)),assigned)]-apply(scores,1,median),
    reference_delta_next=pred$delta.next, reference_best_score=scores[cbind(seq_len(nrow(scores)),assigned)])
}
set.seed(20261114)
common <- readLines(file.path(work,'identity_genes.txt'))
genes <- read.delim(file.path(refdir,'reference_genes.tsv'))$gene_id
rmeta <- read.delim(file.path(refdir,'reference_metadata.tsv'))
cmeta <- read.delim(file.path(refdir,'calibration_metadata.tsv'))
stopifnot(length(common)>5000,all(common %in% genes),!anyDuplicated(common),!any(rmeta$donor_id %in% cmeta$donor_id))
ref <- read_reference('reference'); cal <- read_reference('calibration')
rownames(ref) <- rownames(cal) <- genes
colnames(ref) <- rmeta$reference_barcode; colnames(cal) <- cmeta$reference_barcode
trained <- trainSingleR(normalize(ref[common,,drop=FALSE]), labels=rmeta$broad_type,
  de.method='wilcox',de.n=50,aggr.ref=TRUE,num.threads=4,BPPARAM=MulticoreParam(2))
prediction <- classifySingleR(normalize(cal[common,,drop=FALSE]),trained,num.threads=4)
calibration <- prediction_frame(prediction)
calibration$truth <- cmeta$broad_type[match(calibration$cell_id,cmeta$reference_barcode)]
write.table(calibration,file.path(work,'reference_calibration.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
saveRDS(trained,file.path(work,'reference_trained.rds'))
query <- as(readMM(file.path(work,'identity_query.mtx')),'dgCMatrix')
cells <- readLines(file.path(work,'identity_cells.txt'))
stopifnot(nrow(query)==length(common),ncol(query)==length(cells),length(cells)>0)
rownames(query) <- common; colnames(query) <- cells
prediction <- classifySingleR(normalize(query),trained,num.threads=4)
result <- prediction_frame(prediction)
write.table(result,file.path(work,'identity.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
jsonlite::write_json(list(cells=length(cells),common_nonMT_genes=length(common),
  reference_calibration_accuracy=mean(calibration$truth==calibration$reference_label),
  reference_donors_disjoint=TRUE,domain='adult nuclear reference applied to whole cells; no validated whole-cell accuracy'),
  file.path(work,'IDENTITY_SUCCESS.json'),auto_unbox=TRUE,pretty=TRUE)
