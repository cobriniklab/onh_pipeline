#!/usr/bin/Rscript

copywriter_input_dir = "./output/picard/"
copywriter_output_dir = "./output/copywriter/"
copywriter_capture_regions = "../rb_pipeline/corrected_agilent_regions.bed"

# load required packages
library(CopywriteR)
library(CopyhelpeR)
library(DNAcopy)
library(BiocParallel)


#preCopywriteR

data.folder <- tools::file_path_as_absolute(file.path(copywriter_input_dir))
preCopywriteR(output.folder = file.path(data.folder), bin.size = 20000, ref.genome = "hg19", "chr")


#BiocParallel
bp.param <- MulticoreParam(workers = 4)
bp.param

#CopywriteR
path <- copywriter_input_dir
pattern <- "^.*removed_duplicates.bam$"
copywriter <- function(path, pattern) {
  samples <- list.files(path = path, pattern = pattern, full.names = TRUE)
  controls <- samples
  sample.control <- data.frame(samples, controls)
  CopywriteR(sample.control = sample.control, 
	destination.folder = file.path(copywriter_output_dir), 
	reference.folder = file.path(data.folder, "hg19_20kb_chr"), 
	capture.regions.file = copywriter_capture_regions, 
	bp.param = bp.param) # 
}
copywriter(path, pattern)

plotCNA(destination.folder = copywriter_output_dir)
