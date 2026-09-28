set.seed(42)
rand_cov <- function(n, minlen = 60, maxlen = 3000) {
    lens <- sample(minlen:maxlen, n, replace = TRUE)
    cov <- RleList(lapply(lens, function(l) Rle(rpois(l, 0.3))))
    names(cov) <- paste0("tx", seq_len(n))
    cov
}

testthat::test_that("profile_bin_means matches the per-transcript loop", {
    cov <- rand_cov(200)
    for (nbins in c(50, 100)) {
        orig <- t(sapply(cov, function(x) {
            clos <- nbins * (round(length(x)/nbins, digits = 0) + 1)
            idx <- as.integer(seq(1, length(x), length.out = clos))
            colMeans(matrix(x[idx], ncol = nbins))
        }))
        testthat::expect_identical(profile_bin_means(cov, nbins), orig)
    }
})

testthat::test_that("profile_windows matches the per-transcript loop", {
    cov <- rand_cov(200, minlen = 120)
    idx_fun <- function(l) {
        rnd <- as.integer(l/2)%%3
        mid <- (as.integer(l/2) - rnd)
        c(seq_len(33), (mid - 17):(mid + 15), (l - 32):l)
    }
    orig <- t(sapply(cov, function(x) as.vector(x)[idx_fun(length(x))]))
    testthat::expect_identical(profile_windows(cov, idx_fun), orig)
    # positions beyond the end fall back to the original behaviour (NA)
    past_end <- function(l) c(1:5, l + 1)
    orig <- t(sapply(cov, function(x) as.vector(x)[past_end(length(x))]))
    testthat::expect_identical(profile_windows(cov, past_end), orig)
})

testthat::test_that("coverage_segments behaves like cov[GRanges]", {
    cov <- rand_cov(50)
    gr <- GRanges(sample(names(cov), 80, replace = TRUE), IRanges(1, 1))
    ln <- elementNROWS(cov)[as.character(seqnames(gr))]
    st <- vapply(ln, function(l) sample(seq_len(l - 10), 1), integer(1))
    ranges(gr) <- IRanges(st, width = pmin(ln - st + 1, 200))
    seg <- coverage_segments(cov, gr)
    testthat::expect_identical(segment_names(seg), names(cov[gr]))
    testthat::expect_identical(profile_bin_means(seg, 50), profile_bin_means(cov[gr], 50))
    testthat::expect_identical(profile_windows(seg, function(l) c(1:5, l)), profile_windows(cov[gr],
        function(l) c(1:5, l)))
})

testthat::test_that("extend_terminal_exons matches the per-transcript loop", {
    grl <- GRangesList(a = GRanges("1", IRanges(c(100, 300), width = 50), strand = "+"), b = GRanges("1",
        IRanges(c(900, 700, 500), width = 30), strand = "-"), c = GRanges(), d = GRanges("2", IRanges(1000,
        width = 10), strand = "-"))
    first <- GRangesList(lapply(grl, function(x) {
        if (length(x) == 0) return(x)
        x[1] <- resize(x[1], width = width(x[1]) + 51, fix = "end")
        x
    }))
    last <- GRangesList(lapply(grl, function(x) {
        if (length(x) == 0) return(x)
        x[length(x)] <- resize(x[length(x)], width = width(x[length(x)]) + 51, fix = "start")
        x
    }))
    testthat::expect_equal(as.list(extend_terminal_exons(grl, "first", 51)), as.list(first))
    testthat::expect_equal(as.list(extend_terminal_exons(grl, "last", 51)), as.list(last))
})

testthat::test_that("codon_sums matches aggregate()", {
    gco <- c("AAA", "AAC", "GGG", "TTT")
    pt <- c(AAA = 1, GGG = 2, AAA = 3, TTT = 0, AAA = 5)
    ag <- aggregate(pt, list(names(pt)), sum)
    orig <- rep(0, length(gco))
    names(orig) <- gco
    orig[ag[, 1]] <- ag[, 2]
    testthat::expect_identical(codon_sums(pt, gco), orig)
})

testthat::test_that("readlength padding and counting helpers", {
    d <- data.frame(reads_25 = c(1, 2), reads_27 = c(3, 4), row.names = c("nucl", "chrM"))
    p <- pad_readlength_cols(d, 24:28)
    testthat::expect_identical(colnames(p), paste0("reads_", 24:28))
    testthat::expect_identical(p$reads_26, c(0, 0))
    testthat::expect_identical(p$reads_27, c(3, 4))
    testthat::expect_identical(readlength_union(d, data.frame(reads_30 = 1)), 25:30)
})
