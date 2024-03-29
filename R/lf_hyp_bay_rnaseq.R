#' Summarized Experiment of RNA-seq Counts
#'
#' RNA was harvested from proliferating lung fibroblasts after 72 h of treatment
#' with molidustat or 0.5% oxygen, separately or together with DMSO as the
#' vehicle control. RNA was extracted using Qiagen RNeasy kit and sequenced
#' at BGI.
#'
#' Reads were aligned with `Rsubread` using release 37 of the
#' `GRCh38.primary_assembly.genome.fa.gz` and annotated with the internal *hg38*
#' NCBI RefSeq database. Entrez gene IDs were mapped to HGNC symbols and gene
#' names using `biomaRt`. The sequencing results are available from the NIH short
#' read archive (SRA) and the analysis code is contained in the `data-raw`
#' directory of the package source.
#'
#' Phenotype data can be accessed with `colData(lf_hyp_bay_se)`:
#'
#'  \describe{
#'   \item{id}{sample id}
#'   \item{experiment}{biological replicate samples}
#'   \item{oxygen}{ambient oxygen during culture}
#'   \item{treatment}{
#'   `DMSO` = DMSO vehicle 0.1% \cr
#'   `BAY` = molidustat 10 μM}
#'  }
#'
#' Simple gene annotation information can be accessed with `rowData(lf_hyp_bay_se)`:
#'
#' \describe{
#'   \item{hgnc_symbol}{gene symbol}
#'   \item{description}{gene name}
#' }
"lf_hyp_bay_rnaseq"
