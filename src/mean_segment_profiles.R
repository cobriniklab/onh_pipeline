library(dplyr)

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

chroms <- unique(my_data$chrom)

my_chroms <- lapply(chroms, mean_segment_profile, my_data = my_data)

my_chroms0 <- rbindlist(my_chroms)

gene_lookup <- my_chroms0 %>% 
  filter(!is.na(seg.mean)) %>% 
  filter(seg.mean > 1 | seg.mean < -1) %>% 
  dplyr::mutate(chrom = paste("chr", chrom, sep=""))

gene_lookup.list <- setNames(split(gene_lookup, seq(nrow(gene_lookup))), rownames(gene_lookup))


ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")

gene_sym_from_coord <- function(peak){
  peak <- getBM(
    attributes=c('entrezgene', 'chromosome_name', 'transcript_start', 'transcript_end'), 
    filters = c('chromosome_name', 'start', 'end' ), 
    values = list(gsub("chr", "", as.character(peak$chrom)), peak$start, peak$end), 
    mart = ensembl) %>% 
    mutate(start = peak$start, end = peak$end)
  
}

biomaRt_result <- lapply(gene_lookup.list, gene_sym_from_coord)

biomaRt_result0 <- rbindlist(biomaRt_result)


write.table(my_chroms0, "./results/20170727_tidy_segmentation_table.csv", sep = ",", row.names= F)

test <- mean_segment_profile(my_data, 1)

test <- as.data.frame(test)

p_stchlk <- ggplot(my_chroms0, aes(x = start, xend = end, y = seg.mean, yend = seg.mean)) + geom_segment(size = 0.25) + facet_grid(chrom ~ .)
