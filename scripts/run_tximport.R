args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tximport)
  library('tidyverse')
  library('fs')
  library('rprojroot')
  library(seuratTools)
})

print(stringtiedir)
print(proj_dir)
print(outrds)

proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

txi_features <- seuratTools::load_counts_from_stringtie(proj_dir)

txi_genes <- txi_features$gene
txi_transcripts <- txi_features$transcript


tpm_meta <- seuratTools::load_meta(proj_dir)

feature_seus <- map(list(gene = txi_genes, transcript = txi_transcripts), seu_from_tximport, tpm_meta)

seuratTools::clustering_workflow(proj_dir, feature_seus)

# saveRDS(list(st = st, sg = sg), file = outrds)

sessionInfo()
date()


