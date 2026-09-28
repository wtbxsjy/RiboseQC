testthat::test_that("kmeans3 is reproducible and leaves the caller's random state alone", {
    x <- c(1, 2, 3, 10, 11, 12, 30, 31, 32, 5, 50, 7, 9)
    set.seed(99)
    before <- .Random.seed
    r1 <- kmeans3(x)
    testthat::expect_identical(.Random.seed, before)
    r2 <- kmeans3(x)
    testthat::expect_identical(r1$cluster, r2$cluster)
    # same clusters whatever the caller's seed
    set.seed(5)
    testthat::expect_identical(kmeans3(x)$cluster, r1$cluster)
    # works on the table objects used for the cutoff calculation
    tb <- table(c(rep(1, 10), rep(2, 3), rep(5, 20), rep(9, 1), rep(12, 7)))
    testthat::expect_identical(kmeans3(tb)$cluster, kmeans3(tb)$cluster)
})

testthat::test_that("calc_cutoffs_from_profiles does not depend on the random seed", {
    set.seed(1)
    # same layout as the 5' end profiles: 50 5'UTR, 99 CDS and 50 3'UTR positions
    prof <- matrix(rpois(60 * 199, 1), nrow = 60)
    prof[, 50 + c(13, 16, 19)] <- prof[, 50 + c(13, 16, 19)] + rpois(60 * 3, 15)
    colnames(prof) <- c(paste0("5_UTR_", 1:50), paste0("CDS_", 1:99), paste0("3_UTR_", 1:50))
    prof <- S4Vectors::DataFrame(prof, check.names = FALSE)
    set.seed(10)
    a <- calc_cutoffs_from_profiles(prof, length_max = 29)
    set.seed(20)
    b <- calc_cutoffs_from_profiles(prof, length_max = 29)
    testthat::expect_identical(a, b)
})
