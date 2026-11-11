suppressPackageStartupMessages({library(DropletUtils); library(SingleR); library(scDblFinder)
    library(SingleCellExperiment); library(Matrix); library(BiocParallel)})
root <- path.expand('~/xiaoe-internship-data')
work <- file.path(root,'_processing/all_data_cleaning_20261111/chin_pair')
samples <- read.delim(file.path(work,'samples.tsv'),check.names=FALSE)
trained <- readRDS(file.path(work,'reference_trained.rds'))
common <- readLines(file.path(work,'reference_common_genes.txt'))
normalize <- function(x) {
    x <- as(x,'dgCMatrix'); dn <- dimnames(x)
    x <- x %*% Diagonal(x=10000/pmax(colSums(x),1)); x@x <- log1p(x@x)
    dimnames(x) <- dn; x
}
result <- bplapply(seq_len(nrow(samples)),function(i) {
    sample <- samples$sample_id[[i]]
    directory <- file.path(work,sample)
    if (file.exists(file.path(directory,'R_V2_SUCCESS.txt'))) return(sample)
    set.seed(20261111+i)
    sce <- read10xCounts(file.path(directory,'raw_10x.h5'),type='HDF5',col.names=TRUE)
    counts <- as(counts(sce),'dgCMatrix'); colnames(counts) <- paste0(sample,':',colnames(counts))
    detected <- diff(counts@p); totals <- colSums(counts)
    mt <- 100*colSums(counts[grepl('^MT-',rowData(sce)$Symbol),,drop=FALSE])/pmax(totals,1)
    edpath <- file.path(directory,'emptydrops.tsv.gz')
    if (file.exists(edpath)) {
        eddf <- read.delim(edpath,check.names=FALSE)
        stopifnot(identical(eddf$barcode,colnames(counts)))
        called <- !is.na(eddf$FDR) & eddf$FDR<=.001
        text <- readLines(file.path(directory,'emptydrops_metadata.txt'))
        line <- which(text=='$niters')
        iterations <- as.integer(sub('^\\[1\\] ','',text[line+1]))
    } else {
        ed <- emptyDrops(counts,lower=100,niters=10000,BPPARAM=SerialParam())
        eddf <- as.data.frame(ed); iterations <- 10000L
        called <- !is.na(eddf$FDR) & eddf$FDR<=.001
        limited <- detected>=500 & detected<=5000 & mt<5 & !called & !is.na(eddf$Limited) & eddf$Limited
        if (any(limited)) {
            set.seed(20261111+i)
            ed <- emptyDrops(counts,lower=100,niters=100000,BPPARAM=SerialParam())
            eddf <- as.data.frame(ed); iterations <- 100000L
            called <- !is.na(eddf$FDR) & eddf$FDR<=.001
        }
        con <- gzfile(edpath,'wt')
        write.table(data.frame(barcode=colnames(counts),eddf,emptydrops_called=called),con,
                    sep='\t',row.names=FALSE,quote=FALSE,na='');close(con)
        capture.output(metadata(ed)[c('lower','niters','alpha','retain')],file=file.path(directory,'emptydrops_metadata.txt'))
    }
    # Raw inputs contain many low-information called droplets. Apply the already
    # agreed gene lower bound before estimating a loading-rate-dependent doublet model.
    model <- called & detected>=500
    stopifnot(sum(model)>100)
    set.seed(20261111+i)
    dbl <- scDblFinder(SingleCellExperiment(list(counts=counts[,model,drop=FALSE])),
                        dbr=NULL,dbr.per1k=.008,BPPARAM=SerialParam(),verbose=FALSE)
    prediction <- data.frame(barcode=colnames(dbl),scDblFinder_score=dbl$scDblFinder.score,
                             scDblFinder_class=as.character(dbl$scDblFinder.class))
    write.table(prediction,file.path(directory,'doublets_v2.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    capture.output(metadata(dbl),file=file.path(directory,'scdblfinder_v2_metadata.txt'))
    query <- called & detected>=500 & detected<=5000 & mt<5
    if (!file.exists(file.path(directory,'identity.tsv'))) {
        pred <- classifySingleR(normalize(counts[common,query,drop=FALSE]),trained,num.threads=2)
        scores <- as.matrix(pred$scores); assigned <- match(pred$labels,colnames(scores))
        df <- data.frame(barcode=rownames(pred),reference_label=pred$labels,
                          reference_pruned=pred$pruned.labels,
                          reference_delta=scores[cbind(seq_len(nrow(scores)),assigned)]-apply(scores,1,median),
                          reference_delta_next=pred$delta.next,
                          reference_best_score=scores[cbind(seq_len(nrow(scores)),assigned)])
        write.table(df,file.path(directory,'identity.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
    }
    write.table(data.frame(sample_id=sample,raw_droplets=ncol(counts),called=sum(called),
                           model_nuclei=sum(model),model_min_genes=500,strict_nuclei=sum(query),
                           emptydrops_iterations=iterations,expected_doublet_rate_per1k=.008,
                           expected_loading_rate=.008*sum(model)/1000,
                           emptydrops_limited=sum(eddf$Limited,na.rm=TRUE)),
                file.path(directory,'R_V2_SUMMARY.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    writeLines('Doublet model revised to >=500 genes; identity and cell calling complete',file.path(directory,'R_V2_SUCCESS.txt'))
    cat('Completed v2 ',sample,': model=',sum(model),'; strict=',sum(query),'\n',sep='')
    sample
},BPPARAM=MulticoreParam(workers=4,RNGseed=20261111,stop.on.error=TRUE))
stopifnot(length(result)==13)
capture.output(sessionInfo(),file=file.path(work,'R_v2_sessionInfo.txt'))
writeLines('All 13 captures complete with revised doublet input',file.path(work,'R_V2_SUCCESS.txt'))
