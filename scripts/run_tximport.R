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
print(organism)

organism = c("Mus_musculus" = "mouse", "Homo_sapiens" = "human")[organism]

proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

txi_features <- seuratTools::load_counts_from_stringtie(proj_dir)

tpm_meta <- seuratTools::load_meta(proj_dir)

feature_seus <- map(txi_features, seu_from_tximport, tpm_meta)

feature_seus <- seuratTools::clustering_workflow(proj_dir, feature_seus, organism = organism)

saveRDS(feature_seus, file = outrds)

sessionInfo()
date()


