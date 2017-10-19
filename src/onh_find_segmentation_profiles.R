#!/usr/bin/Rscript

# load required libraries -------------------------------------------------


library(reshape2)
library(tidyverse)
library(data.table)
library(biomaRt)

# add required functions --------------------------------------------------

stepper <- function(fun_chain, x, FUN = print) {
  f_list <- functions(fun_chain)
  for(i in seq_along(f_list)) {
    x <- f_list[[i]](x)
    FUN(x)
  }
  invisible(x)
}

#stepper(my_chain, 0, function(x) {print(x); readline()})

# load data ---------------------------------------------------------------

segmentation_file = list.files("./output/copywriter/CNAprofiles", "^segment.Rdata$", recursive = TRUE, full.names = TRUE)

load(segmentation_file)

# tidy data ---------------------------------------------------------------


format_SCNA_data <- function(dna_copy_object, samples, sample_type = NULL, range.CNA = c(-2, 2),
                          color.palette = colorRampPalette(c("blue", "white", "red"))) {

  
  if(sample_type=="proband"){
    samples <- unique(grep("[[:digit:]]_1", dna_copy_object$output$ID, value = TRUE ))
  } else if (sample_type == "all"){
    samples <- unique(dna_copy_object$output$ID)
  }
  
  samples <- unique(grep("none", samples, value = TRUE))
  
  ## Use all samples by default
  if(missing(samples)) {
    samples <- unique(dna_copy_object$output$ID)
  }
  
  ## Select samples
  dna_copy_object$output <- dna_copy_object$output[dna_copy_object$output$ID %in% samples, ]

  order.samples <- unique(dna_copy_object$output$ID)
  dna_copy_object$output$ID <- factor(dna_copy_object$output$ID, levels = order.samples)
  return(dna_copy_object)

}


SCNA_data <- format_SCNA_data(segment.CNA.object, sample_type="all")
seg_data <- SCNA_data$output
point_data <- gather(SCNA_data$data, sample, value, starts_with("log2"))

chroms <- unique(seg_data$chrom)

lookup_genes <- function(my_data, outfile = NULL) {
  my_data <- my_data %>% 
    rename(sample = ID, start = loc.start, end = loc.end, seg.mean = seg.mean, chrom = chrom)
  
  gene_lookup <- my_data %>% 
    filter(!is.na(seg.mean))
  
  # %>% filter(seg.mean > 0.4 | seg.mean < -0.8) 
  
  gene_lookup.list <- setNames(split(gene_lookup, seq(nrow(gene_lookup))), rownames(gene_lookup))
  
  ensembl  <- useMart(host="grch37.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  
  gene_sym_from_coord <- function(peak){

    peak <- getBM(
      attributes=c('entrezgene', 'hgnc_symbol', 'chromosome_name', 'transcript_start', 'transcript_end'), 
      filters = c('chromosome_name', 'start', 'end' ), 
      values = list(as.character(peak$chrom), peak$start, peak$end), 
      mart = ensembl) %>% 
      mutate(sample = peak$sample, start = peak$start, end = peak$end, seg.mean = peak$seg.mean) %>% 
      filter(!is.na(hgnc_symbol))
    
  }
  
  peak_genes <- lapply(gene_lookup.list, gene_sym_from_coord)
  
  peak_genes <- rbindlist(peak_genes)
  
  # peak_genes <- peak_genes %>% 
  #   dplyr::mutate(sample = gsub(remove_pat, "", sample))
  #   arrange(seg.mean)
  
  write.table(peak_genes, outfile, sep = ",", row.names = FALSE)
  
  return(peak_genes)
  
  
}

segmentation_table <- lookup_genes(seg_data, outfile = "./results/SCNA/SCNA_multisample_peak_genes_all_corrected.csv")


# calculate mean segments for all samples and summarize data --------------

mean_segment_profile <- function(my_data, chromosome){
  my_data <- my_data %>% 
    dplyr::filter(chrom == chromosome) %>% 
    tidyr::gather("loc", "event", 3:4) %>% 
    dplyr::arrange(event) 
  
  out_df <- matrix(nrow=length(my_data$event)/2, ncol=3)
  last_pointer <- min(my_data$event)
  active <- setNames(rownames(my_data[1,]), my_data[1, "ID"])
  start <- my_data$event[[1]]
  for(i in seq_along(my_data$ID)){
    
    pointer <- my_data$event[[i]]
    if((length(active) > 0)&(last_pointer < pointer)){
      
      if(my_data$loc[i] == "loc.start"){
        active <- c(active, setNames(rownames(my_data[i,]), my_data[i,"ID"]))
        
      }
      else{
        end <- my_data$event[[i]]
        set_mean <- mean(my_data$seg.mean[rownames(my_data) %in% active])
        out_df[i/2,] <- c(start, end, set_mean)
        active <- active[!names(active) %in% my_data[i, "ID"]]
        start <- my_data$event[[i]]
      } 
      
    }
  }
  out_df <- as.data.frame(out_df)
  out_df <- out_df %>%
    setNames(c("start", "end", "seg.mean")) %>% 
    mutate(chrom = chromosome)
  return(out_df)
}


mean_chroms <- lapply(chroms, mean_segment_profile, my_data = seg_data)

mean_chroms <- rbindlist(mean_chroms)

write.table(mean_chroms, "./results/SCNA/20170803_mean_segmentation_table.csv", sep = ",", row.names= F)


p_stchlk <- ggplot(my_chroms0, aes(x = start, xend = end, y = seg.mean, yend = seg.mean)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .)


# plot segmentation data --------------------------------------------------

p_seg_stachelek <- ggplot(seg_data, aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean, color = ID)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .) + theme(text = element_text(size=5), legend.position="none")

print(p_seg_stachelek)
ggsave("./doc/stachelek_scna_tumor_segment_plot.png")

p_all <-  ggplot(multi_seg_plot_all, aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean, color = author)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .)
print(p_all)
ggsave("./doc/two_study_compar_scna_tumor_segment_plot.png")

p_line_all <-  ggplot(multi_seg_plot_all, aes(x = loc.start, y = seg.mean, color = author)) + geom_line(size=0.25) + facet_grid(chrom ~ .)
print(p_line_all)
ggsave("./doc/two_study_compar_scna_tumor_segment_plot.png")

