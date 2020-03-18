#!/usr/bin/env Rscript

# load required packages
library(CopywriteR)
library(CopyhelpeR)
library(DNAcopy)
library(BiocParallel)
library('tidyverse')
library('fs')
library('rprojroot')
library('glue')

if (is.null(threads)) threads = 6
if (is.null(copywriter_output_dir)) mypath = "../../output/copywriter"
if (is.null(sample_files))   sample_files <- fs::dir_ls(mypath, glob = "*41-CL_1_recalibrated_10_intronic.bam")
if (is.null(control_files))   control_files <- sample_files 
if (is.null(bin_size)) bin_size = 50000
if (is.null(suffix)) suffix = ""
if (is.null(copywriter_capture_regions)) copywriter_capture_regions = "~/rb_pipeline/corrected_agilent_regions.bed" #must be a .bed file!

#BiocParallel
bp.param <- MulticoreParam(workers = threads)

## ------------------------------------------------------------------------

run_copywriter <- function(bin_size, sample_files, control_files, out_dir, baits_file, ...){
  # browser()
  
  if(length(sample_files) == 0 | !is_file(sample_files)){
    error("please provide sample_files that exist")
  }
  
  dir.create(out_dir)
  
  #preCopywriteR
  
  # data.folder <- tools::file_path_as_absolute(file.path(in_dir))
  preCopywriteR(output.folder = file.path(out_dir), bin.size = bin_size, ref.genome = "hg19", "chr")
  

  
  if(!all(sample_files == control_files)){
    sample_files <- c(sample_files, control_files)
    
    control_files <- c(control_files, control_files)
  }
  

  sample.control <- data.frame(sample = c(sample_files), 
                               control = c(control_files))
  
  CopywriteR(
    sample.control = sample.control,
    destination.folder = file.path(out_dir),
    reference.folder = file.path(out_dir, paste0("hg19_", humanreadable(bin_size), "_chr")),
    bp.param = bp.param,
    capture.regions.file = baits_file,
    ...
  )
  
  #segment and visualize results
  plotCNA(destination.folder = file.path(out_dir))
  
}

## ------------------------------------------------------------------------

# debug(CopywriteR::plotCNA)

run_copywriter(bin_size, sample_files, control_files, copywriter_output_dir, copywriter_capture_regions, keep.intermediary.files = T)

