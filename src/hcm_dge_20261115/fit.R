#!/usr/bin/env Rscript
# Fractions/LODO ranges below describe stability; they are not confidence intervals.
# Use: Rscript fit.R WORK discovery | Rscript fit.R WORK validation [locked_tsv].
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE, digits=16, warn=1)
args <- commandArgs(TRUE)
if (length(args)<2L || !args[2] %in% c("discovery","validation")) stop("Usage: Rscript fit.R WORK discovery|validation [locked_tsv]")
work <- normalizePath(args[1], mustWork=TRUE); phase <- args[2]
pb <- file.path(work,"pb"); out <- file.path(work,"stats")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_out <- function(obj,name) {
  file <- file.path(out,name); tmp <- paste0(file,".",Sys.getpid(),".partial")
  write.table(obj,tmp,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  if(!file.rename(tmp,file)) stop("Could not finalize ",file)
}
read_tsv <- function(path) read.delim(path,check.names=FALSE,row.names=NULL,quote="",comment.char="",na.strings=c("NA",""))
read_pb <- function(modality) {
  path <- file.path(pb,paste0(modality,".tsv.gz"))
  con <- gzfile(path,"rt"); on.exit(close(con))
  d <- read_tsv(con); meta <- read_tsv(file.path(pb,paste0(modality,"_donors.tsv")))
  stopifnot(names(d)[1]=="gene_id",!anyNA(d[[1]]),!anyDuplicated(d[[1]]),
            !anyDuplicated(names(d)),all(c("donor_id","disease") %in% names(meta)),
            !anyNA(meta$donor_id),!anyDuplicated(meta$donor_id),
            setequal(names(d)[-1],as.character(meta$donor_id)))
  genes <- as.character(d[[1]]); x <- as.matrix(d[,-1,drop=FALSE]); rownames(x) <- genes
  stopifnot(is.numeric(x),all(is.finite(x)),all(x>=0),all(x==floor(x)),all(colSums(x)>0))
  m <- meta[match(colnames(x),meta$donor_id),,drop=FALSE]
  m$donor_id <- as.character(m$donor_id)
  m$disease <- factor(as.character(m$disease),levels=c("Control","HCM"))
  stopifnot(!anyNA(m$disease),all(table(m$disease)>0),identical(colnames(x),m$donor_id))
  list(x=x,m=m)
}
design_for <- function(m,adjust) {
  if(adjust) {
    stopifnot(all(c("age","sex") %in% names(m)))
    age <- as.numeric(as.character(m$age)); stopifnot(all(is.finite(age)))
    sex <- tolower(trimws(as.character(m$sex))); sex[sex=="f"] <- "female"; sex[sex=="m"] <- "male"
    stopifnot(all(sex %in% c("female","male"))); m$sex <- factor(sex,levels=c("female","male"))
    m$age_centered <- age-mean(age)
    design <- model.matrix(~age_centered+sex+disease,m)
  } else design <- model.matrix(~disease,m)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design),"diseaseHCM" %in% colnames(design))
  rownames(design) <- m$donor_id
  design
}
keep_primary <- function(obj) {
  y <- DGEList(obj$x)
  rownames(obj$x)[filterByExpr(y,group=obj$m$disease,min.count=10,min.total.count=15,large.n=10,min.prop=.7)]
}
fit_one <- function(obj,genes,adjust=FALSE,skip=integer(),full=TRUE,with_cpm=FALSE) {
  use <- setdiff(seq_len(ncol(obj$x)),skip); m <- obj$m[use,,drop=FALSE]
  genes <- intersect(genes,rownames(obj$x))
  x <- obj$x[genes,use,drop=FALSE]
  x <- x[rowSums(x)>0,,drop=FALSE]; stopifnot(nrow(x)>1,all(colSums(x)>0))
  design <- design_for(m,adjust)
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
  libs <- data.frame(donor_id=m$donor_id,disease=m$disease,source_library_size=colSums(obj$x[,use,drop=FALSE]),
                     analysis_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,
                     effective_library_size=y$samples$lib.size*y$samples$norm.factors)
  list(tab=tab,libs=libs,design=data.frame(donor_id=rownames(design),design,check.names=FALSE))
}
sha256 <- function(file) {
  ret <- system2("sha256sum",shQuote(normalizePath(file,mustWork=TRUE)),stdout=TRUE)
  if(!is.null(attr(ret,"status")) || length(ret)!=1L) stop("Cannot hash ",file)
  substr(ret,1,64)
}
align_table <- function(ids,tab) {
  res <- tab[match(ids,tab$gene_id),,drop=FALSE]; res$gene_id <- ids; rownames(res) <- NULL; res
}
descriptive <- function(obj,ids) {
  # Disease and study/technology are inseparable: equal donor weight, NO hypothesis test.
  kept <- keep_primary(obj)
  y <- calcNormFactors(DGEList(obj$x[kept,,drop=FALSE]),method="TMM")
  log_cp <- log2(cpm(y,log=FALSE,normalized.lib.sizes=TRUE)+0.5)
  c0 <- apply(log_cp[,obj$m$disease=="Control",drop=FALSE],1,median)
  h1 <- apply(log_cp[,obj$m$disease=="HCM",drop=FALSE],1,median)
  res <- data.frame(gene_id=ids,source_feature_present=ids %in% rownames(obj$x),passes_expression_filter=ids %in% kept,
                    median_log2CPM_HCM=unname(h1[ids]),median_log2CPM_Control=unname(c0[ids]),
                    descriptive_log2_cpm_difference=unname((h1-c0)[ids]),stringsAsFactors=FALSE)
  res
}
if(phase=="discovery") {
  raw <- read_pb("scp_raw"); cb <- read_pb("scp_cb")
  stopifnot(identical(rownames(raw$x),rownames(cb$x)),identical(colnames(raw$x),colnames(cb$x)),
            identical(raw$m[,c("donor_id","disease","sex","age")],cb$m[,c("donor_id","disease","sex","age")]),
            all(as.integer(table(raw$m$disease))==c(16L,15L)))
  tested <- keep_primary(raw)
  cat("SCP1303 primary model:",length(tested),"genes across",nrow(raw$m),"donors\n")
  primary <- fit_one(raw,tested,TRUE,with_cpm=TRUE); tab <- primary$tab
  write_out(primary$libs,"scp_raw_libraries.tsv"); write_out(primary$design,"scp_design.tsv")
  sensitivity <- fit_one(cb,tested,TRUE)
  cbtab <- align_table(tab$gene_id,sensitivity$tab)
  tab$cb_logFC <- cbtab$logFC; tab$cb_FDR <- cbtab$FDR
  write_out(sensitivity$libs,"scp_cb_libraries.tsv")
  wk <- max(1L,min(4L,as.integer(Sys.getenv("HCM_LODO_WORKERS","2"))))
  if(.Platform$OS.type=="windows") wk <- 1L
  folds <- parallel::mclapply(seq_len(nrow(raw$m)),function(i) {
    cat("Leave one out:",raw$m$donor_id[i],"\n")
    fit_one(raw,tested,TRUE,skip=i,full=FALSE)[tab$gene_id]
  },mc.cores=wk,mc.preschedule=TRUE)
  if(any(vapply(folds,inherits,logical(1),"try-error"))) stop("A leave-one-donor-out refit failed")
  mat <- do.call(cbind,folds); colnames(mat) <- raw$m$donor_id
  stopifnot(all(is.finite(mat)))
  tab$loo_same_direction_fraction <- rowMeans(sign(mat)==sign(tab$logFC))
  tab$loo_min_abs_logFC <- apply(abs(mat),1,min)
  tab$loo_min_logFC <- apply(mat,1,min); tab$loo_max_logFC <- apply(mat,1,max)
  deltas <- abs(mat-tab$logFC); tab$loo_max_abs_delta <- apply(deltas,1,max)
  tab$loo_max_influence_omitted_donor <- colnames(mat)[max.col(deltas,ties.method="first")]
  loo_path <- file.path(out,"scp_loo_logFC.tsv.gz")
  tmp <- paste0(loo_path,".",Sys.getpid(),".partial")
  con <- gzfile(tmp,"wt")
  write.table(data.frame(gene_id=tab$gene_id,mat,check.names=FALSE),con,sep="\t",row.names=FALSE,quote=FALSE)
  close(con); stopifnot(file.rename(tmp,loo_path))
  tab$hhatl_logFC <- NA_real_; tab$hhatl_FDR <- NA_real_
  reason <- "subtype_matrix_unavailable"
  if(file.exists(file.path(pb,"scp_hhatl.tsv.gz"))) {
    hh <- read_pb("scp_hhatl"); stopifnot("nuclei" %in% names(hh$m))
    elig <- which(hh$m$nuclei>=50); hh$m <- hh$m[elig,,drop=FALSE]; hh$x <- hh$x[,elig,drop=FALSE]
    if(all(table(hh$m$disease)>=5) && !inherits(try(design_for(hh$m,TRUE),silent=TRUE),"try-error")) {
      hfit <- fit_one(hh,tested,TRUE)
      ht <- align_table(tab$gene_id,hfit$tab)
      tab$hhatl_logFC <- ht$logFC; tab$hhatl_FDR <- ht$FDR
      write_out(hfit$libs,"scp_hhatl_libraries.tsv"); reason <- "fitted_on_eligible_donors"
    } else reason <- "insufficient_group_support_or_rank_deficient_design"
  }
  write_out(data.frame(status=reason),"scp_hhatl_status.tsv")
  write_out(tab,"discovery.tsv")
  cat("Discovery complete: ",nrow(tab)," tested; LODO refits ",ncol(mat),"; subtype ",reason,"\n",sep="")
} else {
  lock_file <- if(length(args)>=3) args[3] else file.path(work,"lock","CANDIDATES.tsv")
  seal <- readLines(paste0(lock_file,".sha256"),warn=FALSE); stopifnot(length(seal)==1L)
  stopifnot(tolower(strsplit(trimws(seal),"[[:space:]]+")[[1]][1])==sha256(lock_file))
  lock <- read_tsv(lock_file); idcol <- intersect(c("gene_id","ensembl_id"),names(lock)); stopifnot(length(idcol)>0L)
  ids <- as.character(lock[[idcol[1]]]); stopifnot(!anyNA(ids),!anyDuplicated(ids))
  # Until the immutable list is verified, no held-out counts are opened.
  liu <- read_pb("liu"); stopifnot(all(as.integer(table(liu$m$disease))==c(2L,10L)))
  kept <- keep_primary(liu); liufit <- fit_one(liu,kept,with_cpm=TRUE)
  liutab <- align_table(ids,liufit$tab)
  liutab$source_feature_present <- ids %in% rownames(liu$x)
  liutab$passes_expression_filter <- ids %in% kept
  for(i in which(liu$m$disease=="Control")) {
    # Only effect sizes are reported from remaining-single-control fits.
    eff <- fit_one(liu,kept,skip=i,full=FALSE)
    liutab[[paste0("logFC_omit_",liu$m$donor_id[i])]] <- unname(eff[ids])
  }
  write_out(liutab,"liu.tsv");write_out(liufit$libs,"liu_libraries.tsv")
  for(name in c("chin","whole")) write_out(descriptive(read_pb(name),ids),paste0(name,".tsv"))
  write_out(data.frame(lock_sha256=sha256(lock_file),locked_genes=length(ids),liu_tested_genomewide=length(kept)),"validation_design.tsv")
  cat("Post-lock comparisons complete for ",length(ids)," locked genes\n",sep="")
}
writeLines(capture.output(sessionInfo()),file.path(out,paste0("session_",phase,".txt")))
