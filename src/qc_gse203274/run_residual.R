suppressPackageStartupMessages({library(Matrix);library(scDblFinder);library(SingleCellExperiment);library(BiocParallel)})
work <- path.expand('~/xiaoe-internship-data/_processing/all_data_cleaning_20261111/GSE203274')
summary <- read.delim(file.path(work,'RESIDUAL_REVIEW_SUMMARY.tsv'),check.names=FALSE)
trigger <- tolower(as.character(summary$residual_model_triggered))=='true'
set.seed(20261111)
for(i in which(trigger)) {
    sample <- summary$sample_id[[i]];directory <- file.path(work,'captures',sample)
    x <- as(readMM(file.path(directory,'residual_model_counts.mtx')),'dgCMatrix')
    cells <- readLines(file.path(directory,'residual_model_cells.txt'))
    stopifnot(ncol(x)==length(cells),ncol(x)==summary$model_nuclei[[i]])
    colnames(x) <- cells
    # Published post-QC size is NOT known original loading; caller is supplementary only.
    dbl <- scDblFinder(SingleCellExperiment(list(counts=x)),dbr=NULL,dbr.per1k=.008,
                      BPPARAM=SerialParam(),verbose=FALSE)
    df <- data.frame(cell_id=colnames(dbl),scDblFinder_score=dbl$scDblFinder.score,
                     scDblFinder_class=as.character(dbl$scDblFinder.class))
    write.table(df,file.path(directory,'doublets.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    capture.output(metadata(dbl),file=file.path(directory,'scdblfinder_metadata.txt'))
    cat('Residual model complete:',sample,'; published >=500-gene nuclei=',ncol(x),'\n')
}
writeLines(paste('Residual screen completed for',sum(trigger),'triggered captures; author doublet QC retained'),
           file.path(work,'RESIDUAL_REVIEW_SUCCESS.txt'))
