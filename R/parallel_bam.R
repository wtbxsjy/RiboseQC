#' @include riboseqc.R
NULL

# Parallel processing of one BAM file over genomic windows.
#
# The mapped reads of each chromosome are split into windows of about
# 'target_reads' reads (from the start positions of the reads). Each window is
# read and processed (MAP) independently; a read is assigned to the window
# containing its start, so that reads crossing window boundaries are counted
# once. Window results are then merged with the same REDUCE function used to
# merge the chunks of a sequential reduceByYield() run.

# Genomic windows (GRanges) tiling the chromosomes with reads, each containing
# the start of about 'target_reads' alignments (among those selected by 'flag')
bam_windows <- function(bam_file, target_reads, flag = scanBamFlag()) {
    aln <- scanBam(bam_file, param = ScanBamParam(flag = flag, what = c("rname", "pos")))[[1]]
    seqlens <- seqlengths(seqinfo(BamFile(bam_file)))
    ok <- !is.na(aln$rname) & !is.na(aln$pos)
    by_chr <- split(aln$pos[ok], as.character(aln$rname[ok]))
    wins <- lapply(names(by_chr), function(chr) {
        p <- sort(by_chr[[chr]])
        n_win <- max(1, round(length(p)/target_reads))
        # window starts at read-start quantiles (a pile of reads at one position is not split)
        starts <- unique(c(1, p[round(seq_len(n_win - 1) * length(p)/n_win) + 1]))
        GRanges(chr, IRanges(starts, c(starts[-1] - 1, seqlens[[chr]])))
    })
    unlist(GRangesList(wins), use.names = FALSE)
}

# Number of reads per window used for region-parallel runs: enough windows to
# balance the load over the workers, but never more reads than 'chunk_size'
region_target_reads <- function(bam_file, chunk_size, n_workers) {
    total <- sum(as.numeric(idxstatsBam(bam_file)$mapped))
    max(1, min(chunk_size, ceiling(total/(2 * n_workers))))
}

# Merge a list of chunk results with REDUCE(x = later chunk, y = earlier chunks),
# pairwise, so that each result is merged about log2(n) times
reduce_pairwise <- function(res, REDUCE) {
    while (length(res) > 1) {
        idx <- seq(1, length(res) - 1, by = 2)
        merged <- lapply(idx, function(i) REDUCE(res[[i + 1]], res[[i]]))
        if (length(res)%%2 == 1) {
            merged <- c(merged, res[length(res)])
        }
        res <- merged
    }
    res[[1]]
}

# Run MAP over the reads of a BAM file and merge the results with REDUCE.
# Without BPPARAM this is reduceByYield() over chunks of the file; with a
# BiocParallelParam the file is processed by genomic windows in parallel.
run_bam_pass <- function(bam_file, chunk_size, param, MAP, REDUCE, BPPARAM = NULL, windows = NULL) {
    if (is.null(BPPARAM)) {
        opts <- BamFile(file = bam_file, yieldSize = chunk_size)
        return(reduceByYield(X = opts, YIELD = function(x) readGAlignments(x, param = param),
            MAP = MAP, REDUCE = REDUCE))
    }
    if (is.null(windows)) {
        windows <- bam_windows(bam_file, region_target_reads(bam_file, chunk_size, bpnworkers(BPPARAM)), 
            bamFlag(param, asInteger = TRUE))
    }
    run_window <- function(i) {
        w <- windows[i]
        p <- param
        bamWhich(p) <- w
        x <- readGAlignments(bam_file, param = p)
        x <- x[start(x) >= start(w)]
        # no reads left after removing reads with insertions/deletions (as MAP does)
        if (!any(!grepl("[ID]", cigar(x)))) {
            return(NULL)
        }
        MAP(x)
    }
    # the workers do not use random numbers: keep the caller's RNG state as it was
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (has_seed) {
        seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", seed, envir = globalenv()), add = TRUE)
    }
    res <- bplapply(seq_along(windows), run_window, BPPARAM = BPPARAM)
    res <- res[!vapply(res, is.null, logical(1))]
    if (length(res) == 0) {
        stop(paste("No usable alignments in", bam_file))
    }
    reduce_pairwise(res, REDUCE)
}

# Position tables: counts of reads (or P-sites) per position, kept as plain
# vectors while the chunks/windows of a BAM file are processed and merged
# (GRanges operations have a large per-call overhead), and converted to a
# GRangesList (one element per read length) once the whole file is done.
# group: read length (character); comp: compartment index; sn/st: seqnames and
# strand codes; pos: position; score: count. 'levels' orders the groups.

# Sum the scores of identical entries, sorted by group, comp, seqnames, strand, position
aggregate_positions <- function(group, comp, sn, st, pos, score, levels = NULL) {
    o <- order(as.numeric(group), comp, sn, st, pos)
    new_pos <- c(TRUE, diff(as.numeric(group[o])) != 0 | diff(comp[o]) != 0 | diff(sn[o]) != 0 | diff(st[o]) != 
        0 | diff(pos[o]) != 0)[seq_along(o)]
    first <- o[new_pos]
    score <- rowsum(score[o], cumsum(new_pos), reorder = FALSE)[, 1]
    names(score) <- NULL
    list(group = group[first], comp = comp[first], sn = sn[first], st = st[first], pos = pos[first], score = score, 
        levels = levels)
}

# Merge two position tables (e.g. of two chunks of a BAM file)
merge_position_tables <- function(a, b) {
    aggregate_positions(c(a$group, b$group), c(a$comp, b$comp), c(a$sn, b$sn), c(a$st, b$st), c(a$pos, b$pos), 
        c(a$score, b$score), levels = unique(c(a$levels, b$levels)))
}

# GRangesList with one element per group (read length), in the order of
# 'levels' (default: increasing read length); positions sorted within groups
position_table_to_grl <- function(pt, seqlevs, seqlens, levels = pt$levels) {
    if (is.null(levels)) {
        levels <- as.character(sort(unique(as.numeric(pt$group))))
    }
    strands <- c("+", "-", "*")
    gr <- GRanges(factor(seqlevs[pt$sn], levels = seqlevs), IRanges(pt$pos, width = 1), strand = factor(strands[pt$st], 
        levels = strands), seqlengths = seqlens)
    gr$score <- pt$score
    split(gr, factor(pt$group, levels = levels))
}

# 5' end positions of reads, per read length, as a position table
read_5p_positions <- function(x) {
    plus <- as.character(strand(x)) == "+"
    aggregate_positions(as.character(mcols(x)$len_adj), rep(1L, length(x)), as.integer(seqnames(x)), 
        match(as.character(strand(x)), c("+", "-", "*")), ifelse(plus, start(x), end(x)), rep(1L, length(x)))
}

# P-sites of a chunk of alignments, per read length, for all compartments.
#
# x: alignments without insertions/deletions, soft clips removed from the
# cigars, with mcols 'len_adj' (read length), 'mapq' and 'MD'.
# rl_cutoffs_comp: per compartment, a table with columns read_length and cutoff.
# The P-site of a read is its (cutoff + 1)-th aligned base from the 5' end,
# following splice junctions. Returns GRangesLists (one element per read
# length) of P-site positions with their counts ('score') for all reads
# (P_sites_all), reads with mapq > 50 (P_sites_uniq) and reads with mapq > 50
# and mismatches in the MD tag (P_sites_uniq_mm), as position tables (see
# position_table_to_grl()). Each read length lists the positions of the
# compartments in turn, each sorted.
compute_psites <- function(x, rl_cutoffs_comp, circs, seqlevs) {
    comps <- names(rl_cutoffs_comp)
    sn <- as.character(seqnames(x))
    ct <- rep(NA_real_, length(x))
    comp_idx <- rep(NA_integer_, length(x))
    rls_all <- character(0)
    for (k in seq_along(comps)) {
        resul <- rl_cutoffs_comp[[comps[k]]]
        if (is.null(resul) || length(resul$read_length) == 0) {
            next
        }
        rls_all <- c(rls_all, as.character(as.numeric(resul$read_length)))
        in_comp <- if (comps[k] == "nucl") !sn %in% circs else sn == comps[k]
        m <- match(mcols(x)$len_adj, as.numeric(resul$read_length))
        sel <- in_comp & !is.na(m)
        ct[sel] <- as.numeric(resul$cutoff[m[sel]])
        comp_idx[sel] <- k
    }
    rls_all <- unique(rls_all)
    keep <- which(!is.na(ct))
    x <- x[keep]
    ct <- ct[keep]
    comp_idx <- comp_idx[keep]
    rl <- mcols(x)$len_adj
    plus <- as.character(strand(x)) == "+"
    # position of the P-site in the read, from its leftmost aligned base
    q <- ifelse(plus, ct + 1, rl - ct)
    blocks <- cigarRangesAlongReferenceSpace(cigar(x), pos = start(x), ops = c("M", "=", "X"))
    bl <- unlist(blocks, use.names = FALSE)
    read_id <- rep(seq_along(x), elementNROWS(blocks))
    cum <- unlist(cumsum(width(blocks)), use.names = FALSE)
    before <- cum - width(bl)
    hit <- before < q[read_id] & cum >= q[read_id]
    pos <- rep(NA_real_, length(x))
    pos[read_id[hit]] <- start(bl)[hit] + q[read_id[hit]] - before[hit] - 1
    # offsets outside of the read: continue past its first/last aligned base
    after_end <- q > rl
    pos[after_end] <- end(x)[after_end] + q[after_end] - rl[after_end]
    before_start <- q < 1
    pos[before_start] <- start(x)[before_start] + q[before_start] - 1
    
    uniq <- mcols(x)$mapq > 50
    uniq_mm <- uniq & grepl(x = mcols(x)$MD, pattern = "\\W|.{3,}")
    sn_code <- match(as.character(seqnames(x)), seqlevs)
    st_code <- match(as.character(strand(x)), c("+", "-", "*"))
    count_positions <- function(sel) {
        sel <- which(sel)
        aggregate_positions(as.character(rl[sel]), comp_idx[sel], sn_code[sel], st_code[sel], pos[sel], 
            rep(1L, length(sel)), levels = rls_all)
    }
    list(P_sites_all = count_positions(rep(TRUE, length(x))), P_sites_uniq = count_positions(uniq), 
        P_sites_uniq_mm = count_positions(uniq_mm))
}
