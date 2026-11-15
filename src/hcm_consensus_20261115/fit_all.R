#!/usr/bin/env Rscript
# Full-axis expansion of the already specified donor-count models.
# Liu is discovery evidence in this revision, not a held-out validation cohort.
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE, digits=16, warn=1)

read_tsv <- function(path) read.delim(path,check.names=FALSE,row.names=NULL,
                                     quote="",comment.char="",na.strings=c("NA",""))
read_pb <- function(pb,modality) {
  con <- gzfile(file.path(pb,paste0(modality,".tsv.gz")),"rt")
  on.exit(close(con))
  d <- read_tsv(con)
  meta <- read_tsv(file.path(pb,paste0(modality,"_donors.tsv")))
  stopifnot(names(d)[1]=="gene_id",!anyNA(d[[1]]),!anyDuplicated(d[[1]]),
            !anyDuplicated(names(d)),all(c("donor_id","disease") %in% names(meta)),
            !anyNA(meta$donor_id),!anyDuplicated(meta$donor_id),
            setequal(names(d)[-1],as.character(meta$donor_id)))
  genes <- as.character(d[[1]])
  x <- as.matrix(d[,-1,drop=FALSE]); rownames(x) <- genes
  stopifnot(is.numeric(x),all(is.finite(x)),all(x>=0),all(x==floor(x)),all(colSums(x)>0))
  m <- meta[match(colnames(x),meta$donor_id),,drop=FALSE]
  m$donor_id <- as.character(m$donor_id)
  m$disease <- factor(as.character(m$disease),levels=c("Control","HCM"))
  stopifnot(!anyNA(m$disease),all(table(m$disease)>0),identical(colnames(x),m$donor_id))
  list(x=x,m=m)
}

design_for <- function(m) {
  design <- model.matrix(~disease,m)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design),
            "diseaseHCM" %in% colnames(design))
  rownames(design) <- m$donor_id
  design
}
keep_primary <- function(obj) {
  y <- DGEList(obj$x)
  rownames(obj$x)[filterByExpr(y,group=obj$m$disease,min.count=10,
                             min.total.count=15,large.n=10,min.prop=.7)]
}
fit_one <- function(obj,genes,skip=integer(),full=TRUE,with_cpm=FALSE) {
  # Identical normalization/filtering/fitting order to hcm_dge_20261115/fit.R.
  use <- setdiff(seq_len(ncol(obj$x)),skip); m <- obj$m[use,,drop=FALSE]
  genes <- intersect(genes,rownames(obj$x))
  x <- obj$x[genes,use,drop=FALSE]
  x <- x[rowSums(x)>0,,drop=FALSE]
  stopifnot(nrow(x)>1,all(colSums(x)>0))
  design <- design_for(m)
  y <- calcNormFactors(DGEList(x),method="TMM")
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  q <- glmQLFTest(fit,coef="diseaseHCM")
  tab <- q$table
  stopifnot(identical(rownames(tab),rownames(x)),all(is.finite(tab$logFC)))
  if(!full) return(setNames(tab$logFC,rownames(tab)))
  tab$FDR <- p.adjust(tab$PValue,method="BH")
  tab <- data.frame(gene_id=rownames(tab),tab,row.names=NULL,check.names=FALSE)
  if(with_cpm) {
    norm <- cpm(y,normalized.lib.sizes=TRUE,log=FALSE)
    for(level in c("Control","HCM")) {
      z <- norm[,m$disease==level,drop=FALSE]
      tab[[paste0(level,"_donors_CPM_ge1")]] <- rowSums(z>=1)
      tab[[paste0(level,"_median_CPM")]] <- apply(z,1,median)
    }
  }
  libs <- data.frame(donor_id=m$donor_id,disease=m$disease,
                    source_library_size=colSums(obj$x[,use,drop=FALSE]),
                    analysis_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,
                    effective_library_size=y$samples$lib.size*y$samples$norm.factors)
  list(tab=tab,libs=libs,design=data.frame(donor_id=rownames(design),design,check.names=FALSE))
}
align_table <- function(ids,tab) {
  res <- tab[match(ids,tab$gene_id),,drop=FALSE]
  res$gene_id <- ids; rownames(res) <- NULL
  res
}
descriptive <- function(obj) {
  # Disease and source/technology remain inseparable; no disease P value is emitted.
  ids <- rownames(obj$x); kept <- keep_primary(obj)
  stopifnot(length(kept)>1L)
  y <- calcNormFactors(DGEList(obj$x[kept,,drop=FALSE]),method="TMM")
  log_cp <- log2(cpm(y,log=FALSE,normalized.lib.sizes=TRUE)+0.5)
  c0 <- apply(log_cp[,obj$m$disease=="Control",drop=FALSE],1,median)
  h1 <- apply(log_cp[,obj$m$disease=="HCM",drop=FALSE],1,median)
  tab <- data.frame(gene_id=ids,source_feature_present=TRUE,
                    passes_expression_filter=ids %in% kept,
                    median_log2CPM_HCM=unname(h1[ids]),
                    median_log2CPM_Control=unname(c0[ids]),
                    descriptive_log2_cpm_difference=unname((h1-c0)[ids]),
                    stringsAsFactors=FALSE)
  libs <- data.frame(donor_id=obj$m$donor_id,disease=obj$m$disease,
                    source_library_size=colSums(obj$x),
                    analysis_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,
                    effective_library_size=y$samples$lib.size*y$samples$norm.factors)
  list(tab=tab,libs=libs,kept=kept)
}
write_out <- function(obj,out,name) {
  file <- file.path(out,name)
  stopifnot(!file.exists(file))
  tmp <- paste0(file,".",Sys.getpid(),".partial")
  con <- if(grepl("\\.gz$",name)) gzfile(tmp,"wt") else file(tmp,"wt")
  tryCatch(write.table(obj,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA"),
           finally=close(con))
  stopifnot(file.rename(tmp,file))
}
sha256 <- function(file) {
  ret <- system2("sha256sum",shQuote(normalizePath(file,mustWork=TRUE)),stdout=TRUE)
  if(!is.null(attr(ret,"status")) || length(ret)!=1L) stop("Cannot hash ",file)
  substr(ret,1,64)
}
compare_previous <- function(previous,full) {
  stopifnot(nrow(previous)==200L,!anyNA(previous$gene_id),!anyDuplicated(previous$gene_id),
            all(previous$gene_id %in% full$gene_id),all(names(previous) %in% names(full)))
  current <- align_table(previous$gene_id,full)
  checks <- lapply(names(previous),function(c) {
    a <- previous[[c]]; b <- current[[c]]
    stopifnot(identical(is.na(a),is.na(b)))
    max_abs <- NA_real_
    if(is.numeric(a)) {
      good <- !is.na(a)
      stopifnot(all(is.finite(a[good])),all(is.finite(b[good])))
      delta <- abs(a[good]-b[good]); max_abs <- if(length(delta)) max(delta) else 0
      pass <- all(delta <= 1e-12+3e-9*abs(a[good]))
    } else pass <- identical(as.character(a),as.character(b))
    if(!pass) stop("Prior result changed: ",c)
    data.frame(column=c,rows=length(a),nonmissing=sum(!is.na(a)),
               max_absolute_difference=max_abs,passed=pass)
  })
  do.call(rbind,checks)
}

main <- function(args=commandArgs(TRUE)) {
  if(length(args)!=1L) stop("Usage: Rscript fit_all.R WORK")
  work <- normalizePath(args[1],mustWork=TRUE)
  pb <- file.path(work,"input","pb"); out <- file.path(work,"stats")
  plan <- file.path(work,"code","DESIGN.md")
  previous_file <- file.path(work,"input","previous","liu.tsv")
  stopifnot(file.exists(plan),file.exists(previous_file))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  outputs <- c("liu_all.tsv.gz","chin_all.tsv.gz","whole_all.tsv.gz","FIT_METRICS.tsv",
               "PREVIOUS_LIU_RECONCILIATION.tsv","session_all.txt")
  if(any(file.exists(file.path(out,outputs)))) stop("Use a new workspace; model outputs already exist")
  before <- list()
  for(mod in c("liu","chin","whole")) {
    for(suffix in c(".tsv.gz","_donors.tsv")) {
      f <- file.path(pb,paste0(mod,suffix)); before[[f]] <- sha256(f)
    }
  }
  before[[previous_file]] <- sha256(previous_file)
  before[[plan]] <- sha256(plan)
  liu <- read_pb(pb,"liu")
  stopifnot(all(as.integer(table(liu$m$disease))==c(2L,10L)),
            identical(sort(liu$m$donor_id[liu$m$disease=="Control"]),
                      c("Liu2023_CL0412","Liu2023_CL2")))
  kept <- keep_primary(liu)
  cat("Liu full-source genes:",nrow(liu$x),"; tested after expression filtering:",length(kept),"\n")
  model <- fit_one(liu,kept,with_cpm=TRUE)
  ids <- rownames(liu$x)
  tab <- align_table(ids,model$tab)
  tab$source_feature_present <- TRUE
  tab$passes_expression_filter <- ids %in% kept
  for(i in which(liu$m$disease=="Control")) {
    effect <- fit_one(liu,kept,skip=i,full=FALSE)
    tab[[paste0("logFC_omit_",liu$m$donor_id[i])]] <- unname(effect[ids])
  }
  stopifnot(identical(tab$gene_id,ids),all(is.na(tab$logFC[!tab$passes_expression_filter])),
            all(!is.na(tab$FDR[tab$passes_expression_filter])))
  reconciliation <- compare_previous(read_tsv(previous_file),tab)
  write_out(tab,out,"liu_all.tsv.gz")
  write_out(model$libs,out,"liu_libraries.tsv")
  write_out(model$design,out,"liu_design.tsv")
  write_out(reconciliation,out,"PREVIOUS_LIU_RECONCILIATION.tsv")
  metrics <- data.frame(dataset="liu",source_genes=nrow(liu$x),tested_genes=length(kept),
                        control_donors=sum(liu$m$disease=="Control"),
                        hcm_donors=sum(liu$m$disease=="HCM"),
                        output_rows=nrow(tab),disease_hypothesis_test=TRUE,
                        model="edgeR robust QL; TMM; ~disease",BH_scope="all_expression_filtered_genes")
  for(mod in c("chin","whole")) {
    obj <- read_pb(pb,mod); d <- descriptive(obj)
    write_out(d$tab,out,paste0(mod,"_all.tsv.gz"))
    write_out(d$libs,out,paste0(mod,"_libraries.tsv"))
    metrics <- rbind(metrics,data.frame(dataset=mod,source_genes=nrow(obj$x),tested_genes=length(d$kept),
                                       control_donors=sum(obj$m$disease=="Control"),
                                       hcm_donors=sum(obj$m$disease=="HCM"),output_rows=nrow(d$tab),
                                       disease_hypothesis_test=FALSE,
                                       model="median donor log2(TMM CPM + 0.5) difference",BH_scope="not_applicable"))
  }
  for(f in names(before)) stopifnot(identical(sha256(f),before[[f]]))
  write_out(metrics,out,"FIT_METRICS.tsv")
  write_out(data.frame(path=names(before),sha256=unlist(before,use.names=FALSE),unchanged=TRUE),
            out,"FIT_INPUT_HASHES.tsv")
  writeLines(capture.output(sessionInfo()),file.path(out,"session_all.txt"))
  cat("All-source comparison completed; previous 200 Liu results reconciled.\n")
}
if(sys.nframe()==0L) main()
