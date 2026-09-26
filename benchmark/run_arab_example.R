# Regression run on the package's Arabidopsis example (has circular ChrC/ChrM, so it
# exercises the organelle-compartment code paths that the human chr20 data does not).
# Usage: Rscript run_arab_example.R <pkg_dir> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(args[1])
out_dir <- args[2]
suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
ext <- file.path(pkg_dir, "inst", "ext_data")
annot_dir <- file.path(out_dir, "annot")
fa <- file.path(annot_dir, "simp_arab.fa")
dir.create(annot_dir, showWarnings = FALSE)
if (!file.exists(fa)) {
    writeLines(readLines(gzfile(file.path(ext, "simp_arab.fa.gz"))), fa)
}
annot <- file.path(annot_dir, "simp_arab.gtf.gz_Rannot")
if (!file.exists(annot)) {
    # the example GTF has no "transcript" rows, so use the annotation prebuilt in ex/
    # and point its genome to the local FASTA file
    GTF_annotation <- get(load(file.path(pkg_dir, "ex", "simp_arab.gtf.gz_Rannot")))
    # serialized with old Bioconductor class definitions (e.g. DataFrame -> DFrame)
    old_si <- GTF_annotation$seqinfo
    GTF_annotation <- lapply(GTF_annotation, function(o) if (isS4(o) && !is(o, "Seqinfo")) updateObject(o) else o)
    Rsamtools::indexFa(fa)
    GTF_annotation$genome <- FaFile_Circ(Rsamtools::FaFile(fa), circularRanges = c("ChrC", "ChrM"))
    si <- seqinfo(GTF_annotation$genome)
    GTF_annotation$seqinfo <- si[old_si@seqnames]
    save(GTF_annotation, file = annot)
}
bams <- file.path(ext, c("simp_arab_root.bam", "simp_arab_shoots.bam"))
if (nzchar(Sys.getenv("RIBOSEQC_PROF"))) Rprof(Sys.getenv("RIBOSEQC_PROF"), interval = 0.02, line.profiling = TRUE)
set.seed(1)
RiboseQC_analysis(annotation_file = annot, bam_files = bams,
                  dest_names = file.path(out_dir, c("root", "shoots")),
                  create_report = FALSE, write_tmp_files = TRUE, normalize_cov = FALSE)
