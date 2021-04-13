# rsubread-alignment.R


# build index -------------------------------------------------------------


# ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz

Rsubread::buildindex(
  basename = "data-raw/hs-index/hs-index",
  reference = "/Volumes/Elements/GRCh38.primary_assembly.genome.fa.gz"
)


# align -------------------------------------------------------------------

path <- "/Volumes/Elements/data_raw/rnaseq_lf_hypoxia_molidustat_BGI_2020/Clean/"

reads_1 <-
  list.files(
    path = path,
    pattern = "_1.fq.gz",
    recursive = TRUE,
    full.names = TRUE
  )

reads_2 <-
  list.files(
    path = path,
    pattern = "_2.fq.gz",
    recursive = TRUE,
    full.names = TRUE
  )

output_files <-
  basename(reads_1) %>%
  stringr::str_extract(pattern = ".+(?=\\_1.fq.gz)")

Rsubread::align(
  index = "data-raw/hs-index/hs-index",
  readfile1 = reads_1,
  readfile2 = reads_2,
  output_file = paste0("data-raw/rsubread-bam/", output_files, "_pe.bam"),
  nthreads = 14
)


# counts ------------------------------------------------------------------

bam_files <-
  list.files(
    path = "data-raw/rsubread-bam",
    pattern = "\\.bam$",
    full.names = TRUE
  )

feature_counts <-
  Rsubread::featureCounts(
    files = bam_files,
    annot.inbuilt = "hg38",
    isPairedEnd = TRUE,
    nthreads = 14
  )


# phenotype data ----------------------------------------------------------

pheno_data <-
  tibble::tibble(
    id = sprintf("%02d", 1:16),
    experiment = rep(1:4, each = 4),
    oxygen = rep(c("21%", "0.5%"), 4, each = 2),
    treatment = rep(c("DMSO", "BAY"), 8)
  ) %>%
  dplyr::mutate(
    oxygen = factor(oxygen, levels = c("21%", "0.5%")),
    treatment = factor(treatment, levels = c("DMSO", "BAY"))
  )

usethis::use_data(pheno_data, overwrite = TRUE)


# count data --------------------------------------------------------------

count_data <- feature_counts$counts
colnames(count_data) <-
  stringr::str_extract(colnames(feature_counts$counts), pattern = ".*(?=_pe\\.bam)")

usethis::use_data(count_data, overwrite = TRUE)


# alignment summary -------------------------------------------------------

aligned_reads <-
  feature_counts$stat %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(-Status, names_to = "sample", values_to = "count") %>%
  tidyr::pivot_wider(names_from = Status, values_from = count) %>%
  dplyr::mutate(sample = stringr::str_extract(sample, pattern = ".*(?=_pe\\.bam)")) %>%
  dplyr::rename_with(tolower)

rsubread_alignment_summary <-
  list.files(
    path = "data-raw/rsubread-bam",
    pattern = ".summary",
    full.names = TRUE
  ) %>%
  rlang::set_names(stringr::str_extract(., pattern = "(?<=rsubread-bam/).*(?=_pe\\.bam.summary)")) %>%
  purrr::map_dfr(readr::read_delim, delim = "\t", col_names = FALSE, .id = "sample") %>%
  tidyr::pivot_wider(names_from = X1, values_from = X2) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::left_join(aligned_reads, by = "sample")

usethis::use_data(rsubread_alignment_summary, overwrite = TRUE)
