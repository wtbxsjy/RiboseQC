# Build the RiboseQC annotation (Rannot) for the nf-core/riboseq chr20 test data.
# Usage: Rscript prepare_annotation.R [data_dir]
args <- commandArgs(trailingOnly = TRUE)
data_dir <- normalizePath(if (length(args)) args[1] else "benchmark/data")
suppressPackageStartupMessages(devtools::load_all(normalizePath(file.path(data_dir, "..", "..")), quiet = TRUE))
prepare_annotation_files(annotation_directory = file.path(data_dir, "annot"),
                         genome_seq = file.path(data_dir, "Homo_sapiens.GRCh38.dna.chromosome.20.fa"),
                         gtf_file = file.path(data_dir, "Homo_sapiens.GRCh38.111_chr20.gtf"),
                         scientific_name = "Homo.sapiens", annotation_name = "chr20test",
                         export_bed_tables_TxDb = FALSE, forge_BSgenome = FALSE, create_TxDb = TRUE)
