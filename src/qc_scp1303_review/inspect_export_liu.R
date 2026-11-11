suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(DropletUtils)})
root <- '/home/zhome/xiaoe-internship-data'
out <- file.path(root,'_processing/all_data_cleaning_20261111/liu2023')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
environment <- new.env()
loaded <- load(file.path(root,'validation/liu2023/HCM_snRNAseq_seuratObj.Rdata'),envir=environment)
objects <- loaded[vapply(loaded,function(name) inherits(environment[[name]],'Seurat'),logical(1))]
if(length(objects)!=1L) stop('Expected exactly one Seurat object; inspect loaded names: ',paste(loaded,collapse=','))
object <- environment[[objects[[1]]]]
if(!'RNA' %in% names(object@assays)) stop('RNA count assay unavailable')
counts <- GetAssayData(object,assay='RNA',layer='counts')
metadata <- object[[]]
stopifnot(identical(colnames(counts),rownames(metadata)), !anyDuplicated(colnames(counts)), !anyDuplicated(rownames(counts)))
values <- counts@x
stopifnot(all(is.finite(values)),all(values>=0))
con <- gzfile(file.path(out,'author_metadata.tsv.gz'),'wt')
write.table(data.frame(barcode=rownames(metadata),metadata,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE,na=''); close(con)
summary <- data.frame(field=names(metadata),unique_values=vapply(metadata,function(x)length(unique(x)),integer(1)),
 examples=vapply(metadata,function(x)paste(head(unique(as.character(x)),20),collapse=' | '),character(1)))
write.table(summary,file.path(out,'METADATA_FIELDS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(object_name=objects[[1]],genes=nrow(counts),cells=ncol(counts),
 min_nonzero=min(values),max_count=max(values),fractional_values=sum(abs(values-round(values))>1e-5)),
 file.path(out,'COUNTS_INSPECTION.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
destination <- file.path(out,'liu_RNA_counts.h5')
if(!file.exists(destination)) write10xCounts(destination,counts,gene.id=rownames(counts),gene.symbol=rownames(counts),version='3',type='HDF5')
writeLines('Liu author RNA counts and metadata exported',file.path(out,'AUTHOR_EXPORT_SUCCESS.txt'))
