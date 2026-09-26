# Compare two sets of RiboseQC outputs (*_results_RiboseQC_all, bedgraphs, P_sites_calcs).
# Usage: Rscript compare_results.R <ref_dir> <new_dir>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(GenomicRanges); library(GenomicAlignments)})
ref_dir <- args[1]; new_dir <- args[2]
norm <- function(x) {
    # GRanges: compare as sorted data.frame, ignoring names/seqinfo details
    if (is(x, "GRanges")) {
        x <- sort(x); names(x) <- NULL
        d <- as.data.frame(x); rownames(d) <- NULL
        return(d)
    }
    if (is(x, "GRangesList")) return(lapply(as.list(x), norm))
    if (is(x, "DataFrameList") || is(x, "SimpleList")) return(lapply(as.list(x), norm))
    if (is(x, "DataFrame")) return(norm(as.data.frame(x, optional = TRUE)))
    if (is(x, "RleList")) return(lapply(as.list(x), as.vector))
    if (is.list(x) && !is.data.frame(x)) return(lapply(x, norm))
    x
}
walk <- function(a, b, path = "res_all") {
    a <- norm(a); b <- norm(b)
    if (is.list(a) && !is.data.frame(a) && is.list(b) && !is.data.frame(b)) {
        na <- names(a); nb <- names(b)
        if (!identical(sort(na), sort(nb))) {
            cat("NAMES DIFFER at", path, ":", setdiff(na, nb), "|", setdiff(nb, na), "\n")
        }
        for (n in intersect(na, nb)) walk(a[[n]], b[[n]], paste0(path, "$", n))
        return(invisible())
    }
    r <- all.equal(a, b, tolerance = 0, check.attributes = TRUE)
    if (!isTRUE(r)) {
        cat("DIFF at", path, ":", head(r, 3), "\n")
    } else if (!identical(a, b)) {
        cat("EQUAL VALUES, NOT IDENTICAL (type/attributes) at", path, "\n")
    } else {
        n_ok <<- n_ok + 1
    }
}
n_ok <- 0
for (f in list.files(ref_dir, pattern = "_results_RiboseQC_all$|_results_RiboseQC$|_for_ORFquant$|_junctions$")) {
    cat("==", f, "\n")
    a <- get(load(file.path(ref_dir, f))); b <- get(load(file.path(new_dir, f)))
    walk(a, b)
}
for (f in list.files(ref_dir, pattern = "bedgraph$|_P_sites_calcs$")) {
    same <- identical(readLines(file.path(ref_dir, f)), readLines(file.path(new_dir, f)))
    cat(if (same) "same " else "DIFFERENT ", f, "\n")
}
cat("identical leaves:", n_ok, "\n")
