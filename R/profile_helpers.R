#' @include riboseqc.R
NULL

# Vectorised helpers used to build per-transcript profiles in RiboseQC_analysis.
# They replace per-transcript sapply() loops over RleLists and return exactly the
# same matrices (same values, types and dimnames).

# Concatenated values of an RleList plus the offset of each element in them
.cov_values <- function(cov) {
    if (is(cov, "cov_segments")) {
        return(cov)
    }
    lens <- as.numeric(elementNROWS(cov))
    list(v = as.vector(unlist(cov, use.names = FALSE)), off = cumsum(c(0, lens))[seq_along(lens)],
        len = lens, names = names(cov))
}

# Same content as cov[gr] (an RleList with one element per range of 'gr', named by
# its seqname) for a coverage RleList 'cov', stored as the concatenated values of
# 'cov' plus the offset and length of each range. This avoids extracting one Rle
# per range, which is slow.
coverage_segments <- function(cov, gr) {
    cv <- .cov_values(cov)
    sn <- as.character(seqnames(gr))
    i <- match(sn, names(cov))
    stopifnot(!anyNA(i), all(start(gr) >= 1), all(end(gr) <= cv$len[i]))
    structure(list(v = cv$v, off = cv$off[i] + start(gr) - 1, len = as.numeric(width(gr)), names = sn),
        class = "cov_segments")
}

# Names of the elements of an RleList or of coverage_segments()
segment_names <- function(cov) {
    if (is(cov, "cov_segments")) cov$names else names(cov)
}

# Extend the first or last exon (in list order) of each non-empty transcript by 'by' nt,
# outwards (5' for the first exon, 3' for the last one)
extend_terminal_exons <- function(grl, which = c("first", "last"), by = 51) {
    which <- match.arg(which)
    ex <- unlist(grl, use.names = FALSE)
    part <- PartitioningByEnd(grl)
    idx <- if (which == "first") start(part) else end(part)
    idx <- idx[width(part) > 0]
    ex[idx] <- resize(ex[idx], width = width(ex[idx]) + by, fix = if (which == "first") "end" else "start")
    relist(ex, grl)
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
    res <- matrix(0, nrow = length(L), ncol = nbins, dimnames = list(cv$names, NULL))
    for (cl in unique(clos)) {
        w <- which(clos == cl)
        # positions only depend on the transcript length: compute once per length
        ul <- unique(L[w])
        idx <- lapply(ul, function(l) as.integer(seq(1, l, length.out = cl)))[match(L[w], ul)]
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
        res <- t(sapply(seq_along(cv$len), function(i) {
            cv$v[cv$off[i] + seq_len(cv$len[i])][idx_fun(cv$len[i])]
        }))
        rownames(res) <- cv$names
        return(res)
    }
    matrix(cv$v[unlist(idx, use.names = FALSE) + rep(cv$off, each = k)], ncol = k, byrow = TRUE,
        dimnames = list(cv$names, NULL))
}

# Map the positions of a GRangesList (GRanges with a 'score' column, one element
# per read length) to transcripts in a single call. Returns the mapped hits
# ("all", carrying the score of each hit) and the same hits split by the list
# element (read length) they come from.
map_scored_to_txs <- function(grl, transcripts, seqlevs, seqlens, strand_plus = FALSE) {
    gr <- unlist(grl, use.names = FALSE)
    mp <- mapToTranscripts(gr, transcripts = transcripts)
    mp$score <- gr$score[mp$xHits]
    seqlevels(mp) <- seqlevs
    seqlengths(mp) <- seqlens
    if (strand_plus) {
        strand(mp) <- "+"
    }
    group <- factor(rep(names(grl), elementNROWS(grl))[mp$xHits], levels = unique(names(grl)))
    list(all = mp, by_group = split(mp, group))
}

# Hits belonging to one group ("all" = every hit)
group_hits <- function(mapped, group) {
    if (group == "all") {
        return(mapped$all)
    }
    mapped$by_group[[group]]
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

# Number of reads overlapping each feature, per read length (one column per value
# of 'readlengths'). Same as calling summarizeOverlaps(mode = "Union",
# inter.feature = FALSE, ignore.strand = FALSE) on the reads of each read length
# and binding the counts, but with a single findOverlaps() call.
count_overlaps_by_readlength <- function(features, reads, readlengths) {
    len <- mcols(reads)$len_adj
    ov <- suppressWarnings(findOverlaps(features, reads, ignore.strand = FALSE))
    counts <- table(factor(queryHits(ov), levels = seq_along(features)), factor(len[subjectHits(ov)], 
        levels = readlengths))
    counts <- matrix(as.integer(counts), nrow = length(features), dimnames = list(names(features), 
        paste("reads", readlengths, sep = "_")))
    counts
}

# Contiguous range of read lengths covering the "reads_<length>" columns of two tables
readlength_union <- function(a, b) {
    rls <- as.integer(sub("^reads_", "", c(colnames(a), colnames(b))))
    seq(min(rls), max(rls))
}

# Pad a table with "reads_<length>" columns (data.frame or DataFrame) with zero
# columns so that it covers the read lengths 'rls'
pad_readlength_cols <- function(d, rls) {
    cols <- paste("reads", rls, sep = "_")
    if (identical(colnames(d), cols)) {
        return(d)
    }
    m <- as.matrix(d)
    # zeros of the same type as the counts (integer read counts stay integer)
    out <- matrix(if (is.integer(m)) 0L else 0, nrow = nrow(m), ncol = length(cols), dimnames = list(rownames(m), 
        cols))
    out[, colnames(m)] <- m
    if (is(d, "DataFrame")) {
        return(DataFrame(out, check.names = FALSE))
    }
    as.data.frame(out)
}

# Pad all read-length-dependent tables of a (chunk) result of the BAM statistics pass
pad_readlength_stats <- function(res, rls) {
    for (nm in c("rld", "rld_unq")) {
        res[[nm]] <- pad_readlength_cols(res[[nm]], rls)
    }
    for (nm in c("reads_summary", "reads_summary_unq")) {
        res[[nm]] <- lapply(res[[nm]], pad_readlength_cols, rls = rls)
    }
    res
}

# Read length after removing soft-clipped bases, i.e. qwidth(x) minus the "S" operations
read_length_no_softclip <- function(x) {
    cigarWidthAlongQuerySpace(cigar(x), after.soft.clipping = TRUE)
}

# Counts of all reads and of the reads flagged in 'is_uniq' overlapping each feature.
# Same as two summarizeOverlaps(mode = "Union", inter.feature = FALSE) calls on
# 'reads' and reads[is_uniq], but computing the overlaps only once.
count_overlaps_all_uniq <- function(features, reads, is_uniq, ignore.strand) {
    if (ignore.strand) {
        # as done by summarizeOverlaps()
        if (is(features, "GRangesList")) {
            r <- unlist(features)
            strand(r) <- "*"
            features@unlistData <- r
        } else {
            strand(features) <- "*"
        }
    }
    ov <- findOverlaps(features, reads, ignore.strand = ignore.strand)
    as_counts <- function(cnt) matrix(cnt, ncol = 1, dimnames = list(names(features), "reads"))
    list(all = as_counts(countQueryHits(ov)), uniq = as_counts(countQueryHits(ov[is_uniq[subjectHits(ov)]])))
}
