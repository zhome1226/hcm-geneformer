suppressPackageStartupMessages({
    library(DropletUtils); library(SingleR); library(scDblFinder)
    library(SingleCellExperiment); library(Matrix); library(BiocParallel)
})
root <- path.expand('~/xiaoe-internship-data')
work <- '/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc'
refdir <- file.path(root, 'processed/SCP1303_review_20261111/identity')
samples <- read.delim(file.path(work, 'samples.tsv'), check.names=FALSE)
normalize <- function(x) {
    x <- as(x, 'dgCMatrix')
    dn <- dimnames(x)
    x <- x %*% Diagonal(x=10000/pmax(colSums(x), 1))
    x@x <- log1p(x@x)
    dimnames(x) <- dn
    x
}
read_reference <- function(prefix) {
    con <- gzfile(file.path(refdir, paste0(prefix, '_counts.mtx.gz')), 'rt')
    on.exit(close(con))
    as(readMM(con), 'dgCMatrix')
}
write_prediction <- function(pred, path) {
    scores <- as.matrix(pred$scores)
    assigned <- match(pred$labels, colnames(scores))
    delta <- scores[cbind(seq_len(nrow(scores)), assigned)]-apply(scores, 1, median)
    df <- data.frame(barcode=rownames(pred), reference_label=pred$labels,
                     reference_pruned=pred$pruned.labels, reference_delta=delta,
                     reference_delta_next=pred$delta.next,
                     reference_best_score=scores[cbind(seq_len(nrow(scores)), assigned)])
    write.table(df, path, sep='\t', row.names=FALSE, quote=FALSE, na='')
    df
}
set.seed(20261114)
if (!file.exists(file.path(work, 'REFERENCE_READY.txt'))) {
    first <- read10xCounts(file.path(work,samples$sample_id[[1]],'raw_10x.h5'),
                          type='HDF5', col.names=TRUE)
    genes <- read.delim(file.path(refdir,'reference_genes.tsv'))$gene_id
    common <- genes[genes %in% readLines(file.path(work,'COMMON_GENE_IDS.txt'))]
    stopifnot(length(common)>25000, !anyDuplicated(common))
    ref <- read_reference('reference'); cal <- read_reference('calibration')
    rmeta <- read.delim(file.path(refdir,'reference_metadata.tsv'))
    cmeta <- read.delim(file.path(refdir,'calibration_metadata.tsv'))
    rownames(ref) <- rownames(cal) <- genes
    colnames(ref) <- rmeta$reference_barcode; colnames(cal) <- cmeta$reference_barcode
    trained <- trainSingleR(normalize(ref[common,,drop=FALSE]), labels=rmeta$broad_type,
                            de.method='wilcox', de.n=50, aggr.ref=TRUE,
                            num.threads=8, BPPARAM=MulticoreParam(4))
    pred <- classifySingleR(normalize(cal[common,,drop=FALSE]), trained, num.threads=8)
    table <- write_prediction(pred,file.path(work,'reference_calibration.tsv'))
    table$truth <- cmeta$broad_type[match(table$barcode,cmeta$reference_barcode)]
    write.table(table,file.path(work,'reference_calibration.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
    writeLines(common,file.path(work,'reference_common_genes.txt'))
    saveRDS(trained,file.path(work,'reference_trained.rds'))
    writeLines('Reference trained and held-out donor predictions complete',file.path(work,'REFERENCE_READY.txt'))
    rm(first,ref,cal,pred); gc()
} else {
    trained <- readRDS(file.path(work,'reference_trained.rds'))
    common <- readLines(file.path(work,'reference_common_genes.txt'))
}
cat('Reference ready; processing 6 released capture matrices\n')
results <- bplapply(seq_len(nrow(samples)), function(i) {
    sample <- samples$sample_id[[i]]
    directory <- file.path(work,sample)
    if (file.exists(file.path(directory,'R_SUCCESS.txt'))) return(sample)
    set.seed(20261114+i)
    sce <- read10xCounts(file.path(directory,'raw_10x.h5'),type='HDF5',col.names=TRUE)
    counts <- as(counts(sce),'dgCMatrix')
    colnames(counts) <- paste0(sample,':',colnames(counts))
    detected <- diff(counts@p)
    mito <- grepl('^MT-',rowData(sce)$Symbol)
    totals <- colSums(counts)
    mt <- 100*colSums(counts[mito,,drop=FALSE])/pmax(totals,1)
    # Every released matrix lacks barcodes with positive total UMI <=100.
    # Do not invent an ambient pool or rebrand these as new EmptyDrops calls.
    stopifnot(sum(totals>0 & totals<=100)==0)
    model <- detected>=500
    scored <- sum(model)>=100
    if(scored) {
        dbl <- scDblFinder(SingleCellExperiment(list(counts=counts[,model,drop=FALSE])),
                            dbr=NULL,dbr.per1k=0.008,BPPARAM=SerialParam(),verbose=FALSE)
        dbl_table <- data.frame(barcode=colnames(dbl),scDblFinder_score=dbl$scDblFinder.score,
                                scDblFinder_class=as.character(dbl$scDblFinder.class))
        capture.output(metadata(dbl),file=file.path(directory,'scdblfinder_metadata.txt'))
    } else {
        dbl_table <- data.frame(barcode=colnames(counts)[model],scDblFinder_score=NA_real_,
                                scDblFinder_class='not_scored_insufficient_nuclei')
    }
    write.table(dbl_table,file.path(directory,'doublets_v2.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    query <- detected>=500 & detected<=5000 & mt<5
    if(any(query)) {
        pred <- classifySingleR(normalize(counts[common,query,drop=FALSE]),trained,num.threads=2)
        write_prediction(pred,file.path(directory,'identity.tsv'))
    } else {
        empty <- data.frame(barcode=character(),reference_label=character(),reference_pruned=character(),
                            reference_delta=numeric(),reference_delta_next=numeric(),reference_best_score=numeric())
        write.table(empty,file.path(directory,'identity.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    }
    write.table(data.frame(sample_id=sample,raw_droplets=ncol(counts),released_candidates=ncol(counts),
                           model_nuclei=sum(model),model_min_genes=500,strict_nuclei=sum(query),doublet_scored=scored,
                           expected_doublet_rate_per1k=if(scored) 0.008 else NA_real_,
                           expected_loading_rate=if(scored) .008*sum(model)/1000 else NA_real_,
                           ambient_positive_droplets_UMI_le100=0,
                           cell_calling_status='unverified_released_candidates_no_ambient_pool'),
                file.path(directory,'R_SUMMARY.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    writeLines('Released-candidate doublet screening and reference annotation complete; no EmptyDrops',
               file.path(directory,'R_SUCCESS.txt'))
    cat('Completed ',sample,'; released=',ncol(counts),'; strict=',sum(query),'\n',sep='')
    sample
},BPPARAM=MulticoreParam(workers=4,RNGseed=20261114,stop.on.error=TRUE))
stopifnot(length(results)==6)
capture.output(sessionInfo(),file=file.path(work,'R_sessionInfo.txt'))
writeLines('All 6 sample analyses complete with >=500-gene doublet input; no EmptyDrops',file.path(work,'R_SUCCESS.txt'))
