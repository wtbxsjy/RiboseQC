# RiboseQC analysis benchmark

Scripts used to measure and verify the speed-ups of `RiboseQC_analysis()`.
Every optimisation was checked against the unmodified code: the result
objects (`*_results_RiboseQC_all`, `*_results_RiboseQC`, `*_for_ORFquant`,
`*_junctions`) and the bedgraph / `*_P_sites_calcs` files must be identical.

## Data

* **nf-core/riboseq test data** (human chr20, GRCh38 / Ensembl 111):
  Ribo-seq BAMs `SRX11780887_chr20.bam` (0.60 M reads) and
  `SRX11780888_chr20.bam` (0.57 M reads) from
  [nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/modules/data/genomics/homo_sapiens/riboseq_expression).
* **large**: `SRX11780887_chr20.bam` merged 10 times (6.0 M reads), to
  exercise the BAM-proportional parts and multi-chunk merging.
* **Arabidopsis example** shipped in `inst/ext_data` (has circular ChrC and ChrM,
  i.e. organelle compartments, and many transcripts without UTRs).

## Running

```sh
bash benchmark/get_test_data.sh                  # -> benchmark/data
Rscript benchmark/prepare_annotation.R benchmark/data
# large BAM (optional)
mkdir -p benchmark/data/large
samtools merge -o benchmark/data/large/SRX11780887_x10.bam $(printf 'benchmark/data/SRX11780887_chr20.bam %.0s' {1..10})
samtools index benchmark/data/large/SRX11780887_x10.bam

# run: <package dir> <output dir> [data dir] [Rprof TRUE/FALSE] [chunk_size]
Rscript benchmark/run_benchmark.R . out_new benchmark/data
RIBOSEQC_BAMS=benchmark/data/large/SRX11780887_x10.bam Rscript benchmark/run_benchmark.R . out_large benchmark/data
RIBOSEQC_CORES=2 Rscript benchmark/run_benchmark.R . out_par benchmark/data   # BAMs in parallel
Rscript benchmark/run_arab_example.R . out_arab

# compare with the outputs of another version (e.g. a worktree of the original code)
Rscript benchmark/compare_results.R out_old out_new
```

`run_benchmark.R` calls `set.seed(1)` first: `calc_cutoffs_from_profiles()`
uses `kmeans()` with random starts, so without a fixed seed two runs of the
same code on the same BAM can report different offsets.

The Arabidopsis GTF in `inst/ext_data` has no `transcript` rows, which the
current `prepare_annotation_files()` requires, so `run_arab_example.R` reuses
the annotation prebuilt in `ex/` (updated to current Bioconductor classes and
pointed to the local FASTA). It sets `normalize_cov = FALSE` because the
normalised bedgraph export fails on this tiny example (in the original code too).

## Results

4-core cloud VM, R 4.4.3, Bioconductor 3.20. "Analysis" is the time spent in
`RiboseQC_analysis()`; "wall" also includes starting R and loading the package
(~13 s). Peak memory is that of the largest single process.

| Dataset | Original | Optimised | Speed-up | Outputs |
|---|---|---|---|---|
| chr20, 2 BAMs (1.17 M reads), analysis | 83.6 s | 37.2 s | 2.2x | identical |
| chr20, 2 BAMs, wall / peak memory | 97.6 s / 1465 MB | 50.6 s / 1466 MB | | |
| 6.0 M reads (2 chunks), analysis | 98.4 s | 62.4 s | 1.6x | identical |
| 6.0 M reads, wall / peak memory | 111.8 s / 2611 MB | 75.9 s / 2713 MB | | |
| Arabidopsis example (organelles), wall | 234.2 s | 68.3 s | 3.4x | see below |
| chr20, 4 chunks per BAM (`chunk_size = 150000`) | | | | identical |

Analyzing BAM files in parallel (`BPPARAM = BiocParallel::MulticoreParam(n)`),
wall time: chr20 2 BAMs 50.6 s -> 44.2 s on 2 cores; 6 M + 2 chr20 BAMs
111.2 s -> 79.7 s on 3 cores (bounded by the largest BAM).

Arabidopsis example: the only difference is the "all" 5' profile (and what is
derived from it) of compartments whose transcripts lack UTRs (ChrC, ChrM). The
original code computed it twice and the second pass reused the coverage of the
last read length, so it held the profile of the last read length instead of
the aggregate. Re-introducing that behaviour makes all 54161 compared values
identical to the original output.

Where the time went in the original code (chr20): 5' metagene profiles (~45%,
three `mapToTranscripts()` calls per read length plus per-transcript
`sapply()`), P-site profiles and codon counts (~25%), per-read-length
`summarizeOverlaps()` in the BAM statistics pass, and a full BAM pass only
used to find the read-length range. With 6 M reads, what remains is mostly
reading the BAM file (`readGAlignments()`), overlaps and coverage.
