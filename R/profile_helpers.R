#' @include riboseqc.R
NULL

# Vectorised helpers used to build per-transcript profiles in RiboseQC_analysis.
# They replace per-transcript sapply() loops over RleLists and return exactly the
# same matrices (same values, types and dimnames).

# Concatenated values of an RleList plus the offset of each element in them
.cov_values <- function(cov) {
    lens <- as.numeric(elementNROWS(cov))
    list(v = as.vector(unlist(cov, use.names = FALSE)), off = cumsum(c(0, lens))[seq_along(lens)],
        len = lens)
}

# Equivalent to t(sapply(cov, function(x) {
#   clos <- nbins * (round(length(x)/nbins, digits = 0) + 1)
#   idx <- as.integer(seq(1, length(x), length.out = clos))
#   colMeans(matrix(x[idx], ncol = nbins))
# }))
profile_bin_means <- function(cov, nbins) {
    cv <- .cov_values(cov)
    L <- cv$len
    clos <- nbins * (round(L/nbins, digits = 0) + 1)
    res <- matrix(0, nrow = length(L), ncol = nbins, dimnames = list(names(cov), NULL))
    for (cl in unique(clos)) {
        w <- which(clos == cl)
        idx <- lapply(L[w], function(l) as.integer(seq(1, l, length.out = cl)))
        vals <- cv$v[unlist(idx, use.names = FALSE) + rep(cv$off[w], each = cl)]
        # each column of this matrix holds one bin of one transcript, as in the original matrix()
        means <- colMeans(matrix(vals, nrow = cl/nbins))
        res[w, ] <- matrix(means, ncol = nbins, byrow = TRUE)
    }
    res
}

# Equivalent to t(sapply(cov, function(x) as.vector(x)[idx_fun(length(x))]))
# where idx_fun returns a fixed number of positions for every transcript.
profile_windows <- function(cov, idx_fun) {
    cv <- .cov_values(cov)
    idx <- lapply(cv$len, idx_fun)
    k <- unique(lengths(idx))
    ok <- length(k) == 1 && all(vapply(seq_along(idx), function(i) all(idx[[i]] >= 1 & idx[[i]] <=
        cv$len[i]), logical(1)))
    if (!ok) {
        # out-of-range positions (NA in the original code): fall back to the per-transcript loop
        return(t(sapply(cov, function(x) as.vector(x)[idx_fun(length(x))])))
    }
    matrix(cv$v[unlist(idx, use.names = FALSE) + rep(cv$off, each = k)], ncol = k, byrow = TRUE,
        dimnames = list(names(cov), NULL))
}

# Map the positions of a GRangesList (GRanges with a 'score' column, one element
# per read length) to transcripts in a single call and return the mapped hits,
# carrying the score and the list element (read length) each hit comes from.
map_scored_to_txs <- function(grl, transcripts, seqlevs, seqlens) {
    gr <- unlist(grl, use.names = FALSE)
    mp <- mapToTranscripts(gr, transcripts = transcripts)
    mp$score <- gr$score[mp$xHits]
    mp$group <- rep(names(grl), elementNROWS(grl))[mp$xHits]
    seqlevels(mp) <- seqlevs
    seqlengths(mp) <- seqlens
    mp
}

# Hits belonging to one group ("all" = every hit)
group_hits <- function(mp, group) {
    if (group == "all") {
        return(mp)
    }
    mp[mp$group == group]
}

# Weighted coverage of the hits belonging to one group
group_coverage <- function(mp, group) {
    mp <- group_hits(mp, group)
    coverage(mp, weight = mp$score)
}

# Sum per-transcript values by codon (names(pt)) into a vector over all codons in gco.
# Same result as aggregate(pt, list(names(pt)), sum) followed by filling a zero vector.
codon_sums <- function(pt, gco) {
    sums <- rowsum(unname(pt), names(pt))
    ps_cntt <- rep(0, length(gco))
    names(ps_cntt) <- gco
    ps_cntt[rownames(sums)] <- sums[, 1]
    ps_cntt
}
