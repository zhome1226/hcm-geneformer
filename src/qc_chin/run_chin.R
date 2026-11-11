suppressPackageStartupMessages({
    library(DropletUtils); library(SingleR); library(scDblFinder)
    library(SingleCellExperiment); library(Matrix); library(BiocParallel)
})
root <- path.expand('~/xiaoe-internship-data')
work <- file.path(root, '_processing/all_data_cleaning_20261111/chin_pair')
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
set.seed(20261111)
if (!file.exists(file.path(work, 'REFERENCE_READY.txt'))) {
    first <- read10xCounts(file.path(work,samples$sample_id[[1]],'raw_10x.h5'),
                          type='HDF5', col.names=TRUE)
    genes <- read.delim(file.path(refdir,'reference_genes.tsv'))$gene_id
    common <- genes[genes %in% rownames(first)]
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
cat('Reference ready; processing 13 captures\n')
results <- bplapply(seq_len(nrow(samples)), function(i) {
    sample <- samples$sample_id[[i]]
    directory <- file.path(work,sample)
    if (file.exists(file.path(directory,'R_V2_SUCCESS.txt'))) return(sample)
    set.seed(20261111+i)
    sce <- read10xCounts(file.path(directory,'raw_10x.h5'),type='HDF5',col.names=TRUE)
    counts <- as(counts(sce),'dgCMatrix')
    colnames(counts) <- paste0(sample,':',colnames(counts))
    detected <- diff(counts@p)
    mito <- grepl('^MT-',rowData(sce)$Symbol)
    totals <- colSums(counts)
    mt <- 100*colSums(counts[mito,,drop=FALSE])/pmax(totals,1)
    ed <- emptyDrops(counts,lower=100,niters=10000,BPPARAM=SerialParam())
    eddf <- as.data.frame(ed)
    called <- !is.na(eddf$FDR) & eddf$FDR<=0.001
    unresolved <- detected>=500 & detected<=5000 & mt<5 & !called & !is.na(eddf$Limited) & eddf$Limited
    # Limited nonsignificant strict candidates need greater Monte Carlo resolution.
    iterations <- 10000L
    if (any(unresolved)) {
        set.seed(20261111+i)
        ed <- emptyDrops(counts,lower=100,niters=100000,BPPARAM=SerialParam())
        eddf <- as.data.frame(ed)
        called <- !is.na(eddf$FDR) & eddf$FDR<=0.001
        iterations <- 100000L
    }
    con <- gzfile(file.path(directory,'emptydrops.tsv.gz'),'wt')
    write.table(data.frame(barcode=colnames(counts),eddf,emptydrops_called=called),con,
                sep='\t',row.names=FALSE,quote=FALSE,na=''); close(con)
    capture.output(metadata(ed),file=file.path(directory,'emptydrops_metadata.txt'))
    model <- called & detected>=500
    stopifnot(sum(model)>100)
    dbl <- scDblFinder(SingleCellExperiment(list(counts=counts[,model,drop=FALSE])),
                        dbr=NULL,dbr.per1k=0.008,BPPARAM=SerialParam(),verbose=FALSE)
    dbl_table <- data.frame(barcode=colnames(dbl),scDblFinder_score=dbl$scDblFinder.score,
                            scDblFinder_class=as.character(dbl$scDblFinder.class))
    write.table(dbl_table,file.path(directory,'doublets_v2.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    capture.output(metadata(dbl),file=file.path(directory,'scdblfinder_v2_metadata.txt'))
    query <- called & detected>=500 & detected<=5000 & mt<5
    stopifnot(sum(query)>0)
    pred <- classifySingleR(normalize(counts[common,query,drop=FALSE]),trained,num.threads=2)
    write_prediction(pred,file.path(directory,'identity.tsv'))
    write.table(data.frame(sample_id=sample,raw_droplets=ncol(counts),called=sum(called),
                           model_nuclei=sum(model),model_min_genes=500,strict_nuclei=sum(query),emptydrops_iterations=iterations,
                           expected_doublet_rate_per1k=0.008,expected_loading_rate=.008*sum(model)/1000,
                           emptydrops_limited=sum(eddf$Limited,na.rm=TRUE)),
                file.path(directory,'R_V2_SUMMARY.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    writeLines('EmptyDrops, full-capture doublet model and reference annotation complete',
               file.path(directory,'R_V2_SUCCESS.txt'))
    cat('Completed ',sample,'; called=',sum(called),'; strict=',sum(query),'\n',sep='')
    sample
},BPPARAM=MulticoreParam(workers=4,RNGseed=20261111,stop.on.error=TRUE))
stopifnot(length(results)==13)
capture.output(sessionInfo(),file=file.path(work,'R_sessionInfo.txt'))
writeLines('All 13 sample analyses complete with >=500-gene doublet input',file.path(work,'R_V2_SUCCESS.txt'))
