suppressPackageStartupMessages({
  library(DropletUtils); library(SingleR); library(scDblFinder)
  library(SingleCellExperiment); library(Matrix); library(BiocParallel)
})
root <- path.expand('~/xiaoe-internship-data')
work <- file.path(root,'_processing/all_data_cleaning_20261111/liu2023')
refdir <- file.path(root,'processed/SCP1303_review_20261111/identity')
stopifnot(file.exists(file.path(work,'PREFLIGHT.json')))
samples <- read.delim(file.path(work,'samples.tsv'),check.names=FALSE)
mapping <- read.delim(file.path(work,'reference_gene_mapping.tsv'),check.names=FALSE)
normalize <- function(x) {
  x <- as(x,'dgCMatrix'); dn <- dimnames(x)
  x <- x %*% Diagonal(x=10000/pmax(colSums(x),1)); x@x <- log1p(x@x)
  dimnames(x) <- dn; x
}
read_reference <- function(prefix) {
  con <- gzfile(file.path(refdir,paste0(prefix,'_counts.mtx.gz')),'rt')
  on.exit(close(con)); as(readMM(con),'dgCMatrix')
}
write_prediction <- function(pred,path) {
  scores <- as.matrix(pred$scores)
  assigned <- match(pred$labels,colnames(scores))
  delta <- scores[cbind(seq_len(nrow(scores)),assigned)]-apply(scores,1,median)
  df <- data.frame(barcode=rownames(pred),reference_label=pred$labels,
                   reference_pruned=pred$pruned.labels,reference_delta=delta,
                   reference_delta_next=pred$delta.next,
                   reference_best_score=scores[cbind(seq_len(nrow(scores)),assigned)])
  write.table(df,path,sep='\t',row.names=FALSE,quote=FALSE,na=''); df
}
set.seed(20261111)
common <- mapping$gene_symbol
if (!file.exists(file.path(work,'REFERENCE_READY.txt'))) {
  genes <- read.delim(file.path(refdir,'reference_genes.tsv'))$gene_id
  ref <- read_reference('reference'); cal <- read_reference('calibration')
  rmeta <- read.delim(file.path(refdir,'reference_metadata.tsv'))
  cmeta <- read.delim(file.path(refdir,'calibration_metadata.tsv'))
  stopifnot(!any(rmeta$donor_id %in% cmeta$donor_id))
  rownames(ref) <- rownames(cal) <- genes
  colnames(ref) <- rmeta$reference_barcode; colnames(cal) <- cmeta$reference_barcode
  ref <- ref[mapping$gene_id,,drop=FALSE]; cal <- cal[mapping$gene_id,,drop=FALSE]
  rownames(ref) <- rownames(cal) <- common
  trained <- trainSingleR(normalize(ref),labels=rmeta$broad_type,
                          de.method='wilcox',de.n=50,aggr.ref=TRUE,
                          num.threads=6,BPPARAM=MulticoreParam(3))
  pred <- classifySingleR(normalize(cal),trained,num.threads=6)
  df <- write_prediction(pred,file.path(work,'reference_calibration.tsv'))
  df$truth <- cmeta$broad_type[match(df$barcode,cmeta$reference_barcode)]
  write.table(df,file.path(work,'reference_calibration.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
  writeLines(common,file.path(work,'reference_common_genes.txt'))
  saveRDS(trained,file.path(work,'reference_trained.rds'))
  writeLines('Unique common symbols, fixed reference training donors and held-out donors',file.path(work,'REFERENCE_READY.txt'))
  rm(ref,cal,pred); gc()
} else {
  trained <- readRDS(file.path(work,'reference_trained.rds'))
}
sce <- read10xCounts(file.path(work,'liu_RNA_counts.h5'),type='HDF5',col.names=TRUE)
allcounts <- as(counts(sce),'dgCMatrix')
meta <- read.delim(gzfile(file.path(work,'author_metadata.tsv.gz')),check.names=FALSE)
stopifnot(identical(colnames(allcounts),meta$barcode),nrow(meta)==55122L)
cat('Reference ready; residual-doublet screening and identity on 13 released prefiltered samples\n')
results <- bplapply(seq_len(nrow(samples)),function(i) {
  sample <- samples$sample_id[[i]]; directory <- file.path(work,sample)
  if (file.exists(file.path(directory,'R_SUCCESS.txt'))) return(sample)
  set.seed(20261111+i)
  sel <- meta$orig.ident==sample
  counts <- allcounts[,sel,drop=FALSE]
  colnames(counts) <- paste0(sample,':',colnames(counts))
  detected <- diff(counts@p); total <- colSums(counts)
  mito <- grepl('^MT-',rownames(counts))
  mt <- 100*colSums(counts[mito,,drop=FALSE])/pmax(total,1)
  model <- detected>=500
  # Published Scrublet >.35 plus mixed-lineage filtering is inherited. Only
  # samples with >=10 AND >=1% strict mixed-program flags receive a supplemental
  # screen. Other samples retain author QC without a mechanical second caller.
  trigger <- tolower(as.character(samples$trigger_residual_doublet_screen[[i]]))=='true'
  if (trigger) {
    dbl <- scDblFinder(SingleCellExperiment(list(counts=counts[,model,drop=FALSE])),
                       dbr=NULL,dbr.per1k=.008,BPPARAM=SerialParam(),verbose=FALSE)
    table <- data.frame(barcode=colnames(dbl),scDblFinder_score=dbl$scDblFinder.score,
                        scDblFinder_class=as.character(dbl$scDblFinder.class))
    capture.output(metadata(dbl),file=file.path(directory,'scdblfinder_metadata.txt'))
  } else {
    table <- data.frame(barcode=colnames(counts)[model],scDblFinder_score=NA_real_,
                        scDblFinder_class='not_rerun_author_QC')
  }
  write.table(table,file.path(directory,'doublets.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  strict <- detected>=500 & detected<=5000 & mt<5
  pred <- classifySingleR(normalize(counts[common,strict,drop=FALSE]),trained,num.threads=2)
  write_prediction(pred,file.path(directory,'identity.tsv'))
  write.table(data.frame(sample_id=sample,released_nuclei=ncol(counts),model_nuclei=sum(model),
                         model_min_genes=500,strict_nuclei=sum(strict),
                         residual_screen_triggered=trigger,
                         heuristic_doublet_rate_per1k=if(trigger) .008 else NA_real_,
                         heuristic_rate=if(trigger) .008*sum(model)/1000 else NA_real_,
                         algorithmic_doublet_calls=sum(table$scDblFinder_class=='doublet')),
              file.path(directory,'R_SUMMARY.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  writeLines('Residual-doublet and reference annotation complete',file.path(directory,'R_SUCCESS.txt'))
  cat('Completed ',sample,'; model=',sum(model),'; strict=',sum(strict),'\n',sep='')
  sample
},BPPARAM=MulticoreParam(workers=3,RNGseed=20261111,stop.on.error=TRUE))
stopifnot(length(results)==13L)
capture.output(sessionInfo(),file=file.path(work,'R_sessionInfo.txt'))
writeLines('All 13 released sample analyses complete',file.path(work,'R_SUCCESS.txt'))
