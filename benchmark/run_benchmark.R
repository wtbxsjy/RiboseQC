# Run RiboseQC_analysis on the nf-core/riboseq test BAMs and report timings.
# Usage: Rscript run_benchmark.R <pkg_dir> <out_dir> [data_dir] [profile(TRUE/FALSE)]
args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(args[1])
out_dir <- args[2]
data_dir <- normalizePath(if (length(args) >= 3) args[3] else file.path(pkg_dir, "benchmark", "data"))
do_prof <- length(args) >= 4 && as.logical(args[4])
chunk_size <- if (length(args) >= 5) as.integer(args[5]) else 5000000L
suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
bams <- file.path(data_dir, c("SRX11780887_chr20.bam", "SRX11780888_chr20.bam"))
annot <- file.path(data_dir, "annot", "Homo_sapiens.GRCh38.111_chr20.gtf_Rannot")
dest <- file.path(out_dir, sub(".bam$", "", basename(bams)))
if (do_prof) Rprof(file.path(out_dir, "Rprof.out"), interval = 0.02, line.profiling = TRUE)
set.seed(1)  # calc_cutoffs_from_profiles uses kmeans() with random starts
tm <- system.time(
    RiboseQC_analysis(annotation_file = annot, bam_files = bams, dest_names = dest,
                      create_report = FALSE, write_tmp_files = TRUE,
                      chunk_size = chunk_size)
)
if (do_prof) Rprof(NULL)
print(tm)
writeLines(format(tm[["elapsed"]]), file.path(out_dir, "elapsed_seconds.txt"))
