#!/usr/bin/env Rscript
# Synthetic numerical-compatibility and missing-value test; never reads study outcomes.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
source(args[1])
legacy <- new.env(parent=globalenv())
functions <- c("read_tsv","design_for","keep_primary","fit_one","align_table","descriptive")
for(e in parse(args[2])) {
  if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
     as.character(e[[2]]) %in% functions) eval(e,envir=legacy)
}
stopifnot(all(vapply(functions,exists,logical(1),envir=legacy,inherits=FALSE)))
set.seed(20261115)
donors <- c("Liu2023_CL0412","Liu2023_CL2",sprintf("HCM%02d",1:10))
meta <- data.frame(donor_id=donors,disease=factor(c(rep("Control",2),rep("HCM",10)),
                                                   levels=c("Control","HCM")))
x <- matrix(rnbinom(600*12,mu=100,size=8),600,12,
            dimnames=list(sprintf("ENSG%011d",1:600),donors))
x[551:600,] <- 0L
x[1:40,3:12] <- x[1:40,3:12]*3L
x[41:70,1:2] <- x[41:70,1:2]*3L
x[550,] <- 0L; x[550,1:2] <- 100L
obj <- list(x=x,m=meta)
kept <- keep_primary(obj)
stopifnot(identical(kept,legacy$keep_primary(obj)),length(kept)==550L)
a <- fit_one(obj,kept,with_cpm=TRUE)
b <- legacy$fit_one(obj,kept,with_cpm=TRUE)
stopifnot(isTRUE(all.equal(a$tab,b$tab,tolerance=1e-12)),
          isTRUE(all.equal(a$libs,b$libs,tolerance=1e-12)),
          identical(a$design,b$design))
tab <- align_table(rownames(x),a$tab)
tab$source_feature_present <- TRUE
tab$passes_expression_filter <- rownames(x) %in% kept
for(i in 1:2) {
  u <- fit_one(obj,kept,skip=i,full=FALSE)
  v <- legacy$fit_one(obj,kept,skip=i,full=FALSE)
  stopifnot(isTRUE(all.equal(u,v,tolerance=1e-12)))
  tab[[paste0("logFC_omit_",donors[i])]] <- unname(u[rownames(x)])
}
stopifnot(all(is.na(tab$FDR[551:600])),all(is.na(tab$logFC[551:600])),
          all(tab$source_feature_present),identical(tab$gene_id,rownames(x)))
prev <- tab[c(200:1),]
check <- compare_previous(prev,tab)
stopifnot(all(check$passed))
bad <- prev; bad$logFC[1] <- bad$logFC[1]+.1
stopifnot(inherits(try(compare_previous(bad,tab),silent=TRUE),"try-error"))
d <- descriptive(obj)
expected <- legacy$descriptive(obj,rownames(x))
stopifnot(isTRUE(all.equal(d$tab,expected,tolerance=1e-12)),
          all(is.na(d$tab$descriptive_log2_cpm_difference[551:600])),
          !any(c("PValue","FDR","logFC") %in% names(d$tab)))
out <- normalizePath(args[3],mustWork=TRUE)
write_out(tab,out,"synthetic_all.tsv.gz")
roundtrip <- read_tsv(gzfile(file.path(out,"synthetic_all.tsv.gz"),"rt"))
stopifnot(isTRUE(all.equal(tab,roundtrip,tolerance=1e-12,check.attributes=FALSE)))
write_out(check,out,"synthetic_reconciliation.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"session_test.txt"))
writeLines("passed: legacy numerical match; both control omissions; full axes and 50 filtered NA rows; reordered 200-row reconciliation; tamper rejection; descriptive-only estimates; gzip roundtrip",
           file.path(out,"TEST_RESULT.txt"))
cat(readLines(file.path(out,"TEST_RESULT.txt")),"\n")
