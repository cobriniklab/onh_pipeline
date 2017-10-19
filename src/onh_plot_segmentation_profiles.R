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


SingleSegPlot <- function(dna_copy_object, samples, sample_type = NULL, range.CNA = c(-2, 2),
                          color.palette = colorRampPalette(c("blue", "white", "red"))) {
  browser()
  
  if(sample_type=="proband"){
    samples <- unique(grep("[[:digit:]]_1", dna_copy_object$output$ID, value = TRUE ))
  } else if (sample_type == "all"){
    samples <- unique(grep)
  }
  
  samples <- unique(grep("none", samples, value = TRUE))
  
  ## Use all samples by default
  if(missing(samples)) {
    samples <- unique(dna_copy_object$output$ID)
  }
  
  ## Select samples
  dna_copy_object$output <- dna_copy_object$output[dna_copy_object$output$ID %in% samples, ]
  ## Cap range
  #dna_copy_object$output$seg.mean[dna_copy_object$output$seg.mean < range.CNA[1]] <- range.CNA[1]
  #dna_copy_object$output$seg.mean[dna_copy_object$output$seg.mean > range.CNA[2]] <- range.CNA[2]

  ## Reshape data.frame according to sample name
  # reshape2::dcast -> to wider data.frame; right-hand of tilde needs to be a 'measured variable' (i.e., needs to go into columns)
  # reshape2::melt -> to narrower data.frame
  # Names are changed by dcast according to level order -> change order levels (!)
  order.samples <- unique(dna_copy_object$output$ID)
  dna_copy_object$output$ID <- factor(dna_copy_object$output$ID, levels = order.samples)
  return(dna_copy_object)
  #dna_copy_object.cast <- reshape2::dcast(data = dna_copy_object$output, formula = ID + chrom + loc.start + loc.end + num.mark + seg.mean ~ ID
  #                                        , value.var = "num.mark")

  #dna_copy_object.cast[is.na(dna_copy_object.cast)] <- 0

  ## Calculate number of samples
  sample.number <- ncol(dna_copy_object.cast) - 6

  ## Collapse data
  #dna_copy_object.cast.aggregate <- aggregate(dna_copy_object.cast[, c("chrom", "num.mark")], by = list(dna_copy_object.cast$chrom), FUN = sum)

  ## Calculate scaling factors
  #range.factor <- range.CNA[2] - range.CNA[1]

  ## Create overviewPlot
  #my_data_2 <-  as.data.frame(dna_copy_object.cast) %>%
   # dplyr::rename_(point_density = names(.)[7])
}


my_data <- SingleSegPlot(segment.CNA.object, sample_type="proband")
seg_data <- my_data$output
point_data <- gather(my_data$data, sample, value, starts_with("log2"))


p_seg_stachelek <- ggplot(seg_data, aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean, color = ID)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .) + theme(text = element_text(size=5), legend.position="none")

print(p_seg_stachelek)
ggsave("./doc/stachelek_scna_tumor_segment_plot.png")


p_all <-  ggplot(multi_seg_plot_all, aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean, color = author)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .)
print(p_all)
ggsave("./doc/two_study_compar_scna_tumor_segment_plot.png")

p_line_all <-  ggplot(multi_seg_plot_all, aes(x = loc.start, y = seg.mean, color = author)) + geom_line(size=0.25) + facet_grid(chrom ~ .)
print(p_line_all)
ggsave("./doc/two_study_compar_scna_tumor_segment_plot.png")

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

chroms <- unique(seg_data$chrom)

mean_chroms <- lapply(chroms, mean_segment_profile, my_data = seg_data)

mean_chroms <- rbindlist(mean_chroms)

write.table(mean_chroms, "./results/SCNA/20170803_mean_segmentation_table.csv", sep = ",", row.names= F)

lookup_genes <- function(my_data) {
  my_data <- my_data %>% 
    rename(sample = ID, start = loc.start, end = loc.end, seg.mean = seg.mean, chrom = chrom)
  
  gene_lookup <- my_data %>% 
    filter(!is.na(seg.mean)) 
  # %>%   filter(seg.mean > 0.4 | seg.mean < -0.8) 
  
  gene_lookup.list <- setNames(split(gene_lookup, seq(nrow(gene_lookup))), rownames(gene_lookup))
  
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  
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
  
  peak_genes <- peak_genes %>% 
    arrange(seg.mean)
  
  return(peak_genes)
  
  write.table(peak_genes, "./results/SCNA/SCNA_multisample_peak_genes_probands.csv", sep = ",", row.names = FALSE)
  
}

test0 <- lookup_genes(seg_data)

test <- mean_segment_profile(my_data, 1)

test <- as.data.frame(test)

p_stchlk <- ggplot(my_chroms0, aes(x = start, xend = end, y = seg.mean, yend = seg.mean)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .)



