suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix); library(DropletUtils)})
root <- path.expand('~/xiaoe-internship-data')
work <- file.path(root,'_processing/all_data_cleaning_20261111/gse138262')
src <- file.path(root,'validation/hcm_GSE138262/author_corrected_zenodo18206554')
results <- list()
for (kind in c('ROOIJonly','ALL')) {
    paths <- list.files(src,pattern=paste0('_',kind,'[.]'),full.names=TRUE)
    stopifnot(length(paths)==1L)
    message('Reading ',basename(paths))
    object <- readRDS(paths)
    stopifnot(inherits(object,'Seurat'),'RNA' %in% names(object@assays))
    counts <- as(GetAssayData(object,assay='RNA',layer='counts'),'dgCMatrix')
    metadata <- object[[]]
    stopifnot(identical(colnames(counts),rownames(metadata)),!anyDuplicated(colnames(counts)),!anyDuplicated(rownames(counts)))
    stopifnot(all(is.finite(counts@x)),all(counts@x>=0))
    directory <- file.path(work,paste0('author_',kind))
    dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    con <- gzfile(file.path(directory,'metadata.tsv.gz'),'wt')
    write.table(data.frame(barcode=rownames(metadata),metadata,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE,na='')
    close(con)
    fields <- data.frame(field=names(metadata),class=vapply(metadata,function(x)paste(class(x),collapse=','),character(1)),
                         unique_values=vapply(metadata,function(x)length(unique(x)),integer(1)),
                         examples=vapply(metadata,function(x)paste(head(unique(as.character(x)),20),collapse=' | '),character(1)))
    write.table(fields,file.path(directory,'METADATA_FIELDS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
    for (assay_name in names(object@assays)) {
        assay_counts <- tryCatch(as(GetAssayData(object,assay=assay_name,layer='counts'),'dgCMatrix'),error=function(e)NULL)
        if(is.null(assay_counts) || nrow(assay_counts)==0L || ncol(assay_counts)==0L) next
        destination <- file.path(directory,paste0(assay_name,'_counts.h5'))
        if(file.exists(destination)) stop('Author export already exists; inspect it before rerunning')
        write10xCounts(destination,assay_counts,gene.id=rownames(assay_counts),gene.symbol=rownames(assay_counts),version='3',type='HDF5')
        write.table(data.frame(assay=assay_name,genes=nrow(assay_counts),cells=ncol(assay_counts),
                              nnz=length(assay_counts@x),fractional_values=sum(abs(assay_counts@x-round(assay_counts@x))>1e-5)),
                    file.path(directory,paste0(assay_name,'_inspection.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
    }
    for (reduction in names(object@reductions)) {
        coords <- Embeddings(object[[reduction]])
        con <- gzfile(file.path(directory,paste0('author_',reduction,'_coordinates.tsv.gz')),'wt')
        write.table(data.frame(barcode=rownames(coords),coords,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE,na='')
        close(con)
    }
    results[[kind]] <- data.frame(kind=kind,source=paths,class=paste(class(object),collapse=','),
                                  genes=nrow(counts),cells=ncol(counts),nnz=length(counts@x),
                                  fractional_values=sum(abs(counts@x-round(counts@x))>1e-5),
                                  gene_examples=paste(head(rownames(counts),6),collapse=' | '),
                                  cell_examples=paste(head(colnames(counts),6),collapse=' | '),
                                  assays=paste(names(object@assays),collapse=' | '))
    rm(object,counts,metadata);gc()
}
write.table(do.call(rbind,results),file.path(work,'AUTHOR_OBJECT_INSPECTION.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines('R RNA counts and complete author metadata exported',file.path(work,'AUTHOR_EXPORT_SUCCESS.txt'))
