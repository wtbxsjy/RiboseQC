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
RIBOSEQC_REGION_CORES=4 Rscript benchmark/run_benchmark.R . out_reg benchmark/data   # each BAM over genomic windows
Rscript benchmark/run_arab_example.R . out_arab   # also takes RIBOSEQC_REGION_CORES

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

### Sequential

| Dataset | Original | Now | Speed-up |
|---|---|---|---|
| chr20, 2 BAMs (1.17 M reads), analysis | 83.6 s | 35.7 s | 2.3x |
| chr20, 2 BAMs, wall / peak memory | 97.6 s / 1465 MB | 48.6 s / 1482 MB | |
| 6.0 M reads (2 chunks), analysis | 98.4 s | 55.3 s | 1.8x |
| 6.0 M reads, wall / peak memory | 111.8 s / 2611 MB | 68.1 s / 2722 MB | |
| Arabidopsis example (organelles), wall | 234.2 s | 68.3 s | 3.4x |

### Parallel

* `BPPARAM = BiocParallel::MulticoreParam(n)`: BAM files in parallel.
  chr20, 2 BAMs on 2 cores: analysis 35.7 s -> 22.0 s.
* `region_BPPARAM = BiocParallel::MulticoreParam(n)`: each BAM file in
  parallel over genomic windows (`R/parallel_bam.R`). The chromosomes are cut
  into windows holding about `total reads / (2 n)` read starts (at most
  `chunk_size`), each read belongs to the window containing its start, and
  window results are merged like the chunks of a sequential run.
  6.0 M reads: analysis 55.3 s -> 53.5 s on 2 cores, 41.2 s on 4 cores, and
  peak memory 2722 MB -> 1598 MB (each process holds one window). No gain on
  the small chr20 BAMs (0.6 M reads), where the BAM passes take only a few
  seconds.

  On this VM, reading BAM files scales poorly with processes
  (`readGAlignments()` on 4 regions: 8.7 s sequentially, 6.4 s on 4 workers),
  which bounds the gain here; the parallel efficiency of plain R code is ~75%.

### Outputs

Sequential runs with one chunk, with 4 chunks (`chunk_size = 150000`) and
region-parallel runs on 2 and 4 cores give identical results (chr20, 6 M reads,
Arabidopsis example). Compared with the original code, the differences are all
fixes:

* `P_sites_uniq_mm` (and `for_ORFquant$P_sites_uniq_mm`): the original took
  the P-sites of the wrong reads on the plus strand (logical index built on the
  reads, applied to P-sites in another order) and whole reads instead of
  P-sites on the minus strand, so the result also depended on the chunking.
  P-sites are now computed per read for all read lengths at once
  (`compute_psites()`); `P_sites_all` and `P_sites_uniq` are identical to the
  original on all four test BAMs.
* `sequence_analysis$len_adj` held values like 28.28 instead of 28 (from the
  names of `unlist()`).
* Merging chunks no longer goes through `coverage()`, which joined adjacent
  positions with equal counts into wider ranges (dropped by
  `mapToTranscripts()` across exon junctions): 5' end and P-site positions are
  counted per position whatever the chunking.
* `reads_summary` counts are always integer (they were double when a read
  length was missing from a chunk), and the minus-strand coverage in
  `P_sites_stats` is always named `coverage_*_minus`.
* Arabidopsis example: the "all" 5' profile of compartments whose transcripts
  lack UTRs (ChrC, ChrM). The original code computed it twice and the second
  pass reused the coverage of the last read length, so it held the profile of
  the last read length instead of the aggregate. Re-introducing that behaviour
  made all 54161 compared values identical to the original output.

Where the time went in the original code (chr20): 5' metagene profiles (~45%,
three `mapToTranscripts()` calls per read length plus per-transcript
`sapply()`), P-site profiles and codon counts (~25%), per-read-length
`summarizeOverlaps()` in the BAM statistics pass, and a full BAM pass only
used to find the read-length range. With 6 M reads, what remains is mostly
reading the BAM file (`readGAlignments()`), overlaps and coverage.
