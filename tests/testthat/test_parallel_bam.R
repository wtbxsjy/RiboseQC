bam <- system.file(package = "RiboseQC", "ext_data", "simp_arab_root.bam")
if (!nzchar(bam)) bam <- file.path("..", "..", "inst", "ext_data", "simp_arab_root.bam")

testthat::test_that("bam_windows tiles every chromosome with reads, without gaps", {
    w <- bam_windows(bam, 2000)
    lens <- seqlengths(seqinfo(Rsamtools::BamFile(bam)))
    for (chr in unique(as.character(seqnames(w)))) {
        wc <- w[seqnames(w) == chr]
        testthat::expect_identical(start(wc)[1], 1L)
        testthat::expect_identical(end(wc)[length(wc)], as.integer(lens[[chr]]))
        testthat::expect_true(all(start(wc)[-1] == end(wc)[-length(wc)] + 1L))
    }
    testthat::expect_gt(length(w), length(unique(seqnames(w))))
})

testthat::test_that("merge_position_tables sums counts of identical positions", {
    a <- aggregate_positions(c("28", "28", "29"), c(1L, 1L, 1L), c(1L, 1L, 2L), c(1L, 1L, 2L), c(10, 10, 5), 
        c(1L, 2L, 1L))
    testthat::expect_identical(a$score, c(3L, 1L))
    b <- aggregate_positions(c("28", "30"), c(1L, 1L), c(1L, 1L), c(1L, 1L), c(10, 7), c(4L, 1L))
    m <- merge_position_tables(a, b)
    testthat::expect_identical(m$group, c("28", "29", "30"))
    testthat::expect_identical(m$score, c(7L, 1L, 1L))
    grl <- position_table_to_grl(m, c("1", "2"), c(`1` = 100, `2` = 100))
    testthat::expect_identical(names(grl), c("28", "29", "30"))
    testthat::expect_identical(sum(unlist(grl)$score), 9L)
})

testthat::test_that("run_bam_pass gives the same result sequentially and by genomic windows", {
    param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), what = "mapq")
    MAP <- function(x) list(n = table(factor(as.character(seqnames(x)), levels = seqlevels(x))), pos = read_5p_positions({
        mcols(x)$len_adj <- qwidth(x)
        x
    }))
    REDUCE <- function(x, y) list(n = x$n + y$n, pos = merge_position_tables(x$pos, y$pos))
    seq_res <- run_bam_pass(bam, 5000, param, MAP, REDUCE)
    win_res <- run_bam_pass(bam, 5000, param, MAP, REDUCE, BPPARAM = BiocParallel::SerialParam(), windows = bam_windows(bam, 
        3000))
    testthat::expect_identical(win_res$n, seq_res$n)
    testthat::expect_identical(win_res$pos[c("group", "sn", "st", "pos", "score")], seq_res$pos[c("group", 
        "sn", "st", "pos", "score")])
})
