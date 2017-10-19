#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(splitstackshape)

##Read files named 41-T_1_5%.hg19_multianno.txt, etc.
laura_li_proband_filenames <- list.files(path="./output/annovar/", pattern="[[:digit:]]_.*ONH.*.csv", full.names = TRUE)
laura_li_parental_filenames <- list.files(path="./output/annovar/", pattern="[[:digit:]]-.*ONH.*.csv" , full.names = TRUE)
haplo_proband_filenames <- list.files(path="./output/annovar/", pattern="[[:digit:]]_1_.*haplo.*.csv", full.names = TRUE)
haplo_parental_filenames <- list.files(path="./output/annovar/", pattern="[[:digit:]]-.*_1_haplo.*.csv", full.names = TRUE)



##Create list of data frame names without the ".txt" part 
laura_li_proband_names <-substr(laura_li_proband_filenames,19,21) #format for laura_li samples
laura_li_parental_names <-substr(laura_li_parental_filenames,19,23) #format for laura_li samples
haplo_proband_names <-substr(haplo_proband_filenames,19,29) #format for haplocaller samples
haplo_parental_names <-substr(haplo_parental_filenames,19,33) #format for haplocaller samples


# # load laura li parentals ---------------------------------------------------
# laura_li_parental_list <- lapply(laura_li_parental_filenames, function(x)read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE))
# names(laura_li_parental_list) <- laura_li_parental_names
# sample_annotation = c("GT", "DP", "EC", "SGCONFS", "SGCOUNTREF_F", "SGCOUNTREF_R", "SGCOUNTALT_F", "SGCOUNTALT_R", "SGBRC", "SGBRP")
# tidy_laura_li_parental <- laura_li_parental_list%>%
#   rbindlist(use.names=TRUE, fill=TRUE, idcol="sample") %>%
#   select(sample, Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, gnomAD_exome_ALL, SIFT_score, Otherinfo) %>%
#   separate(Otherinfo, c("a", "b", "c", "d", "e", "f", "g", "h", "i", "filter_status", "k", "annotation", "annotation_results"), "\t") %>%
#   select(-a, -b, -c, -d, -e, -f, -g, -h, -i, -k, -annotation) %>%
#   separate(annotation_results, sample_annotation, ":")


# load laura li probands ----------------------------------------------------
laura_li_proband_list <- lapply(laura_li_proband_filenames, function(x)read.table(x, header = TRUE, sep = ",", stringsAsFactors = FALSE, fill=TRUE))
names(laura_li_proband_list) <- laura_li_proband_names
sample_annotation = c("GT", "DP", "EC", "SGCONFS", "SGCOUNTREF_F", "SGCOUNTREF_R", "SGCOUNTALT_F", "SGCOUNTALT_R", "SGBRC", "SGBRP")
tidy_laura_li_proband <- laura_li_proband_list%>%
  rbindlist(use.names=TRUE, fill=TRUE, idcol="sample") %>%
  select(sample, Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, gnomAD_exome_ALL, SIFT_score, Otherinfo) %>%
  separate(Otherinfo, c("a", "b", "c", "d", "e", "f", "g", "h", "i", "filter_status", "k", "annotation", "annotation_results"), "\t") %>%
  select(-a, -b, -c, -d, -e, -f, -g, -h, -i, -k, -annotation) %>%
  separate(annotation_results, sample_annotation, ":")


# load haplotypecaller parentals ------------------------------------------
# haplotypecaller_parental_list <- lapply(haplotypecaller_parental_filenames, function(x)read.table(x, header = TRUE, sep = ",", stringsAsFactors = FALSE, fill=TRUE))
# names(haplotypecaller_parental_list) <- haplotypecaller_parental_names
# sample_annotation = c("GT", "DP", "EC", "SGCONFS", "SGCOUNTREF_F", "SGCOUNTREF_R", "SGCOUNTALT_F", "SGCOUNTALT_R", "SGBRC", "SGBRP")
# tidy_haplo_parental <- haplotypecaller_parental_list%>%
#   rbindlist(use.names=TRUE, fill=TRUE, idcol="sample") %>%
#   select(sample, Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, gnomAD_exome_ALL, SIFT_score, Otherinfo) %>%
#   separate(Otherinfo, c("a", "b", "c", "d", "e", "f", "g", "h", "i", "filter_status", "k", "annotation", "annotation_results"), "\t") %>%
#   select(-a, -b, -c, -d, -e, -f, -g, -h, -i, -k, -annotation) %>%
#   separate(annotation_results, sample_annotation, ":")

# load haplotypecaller probands -------------------------------------------
haplo_proband_list <- lapply(haplo_proband_filenames, function(x)read.table(x, header = TRUE, sep = ",", stringsAsFactors = FALSE, fill=TRUE))
names(haplo_proband_list) <- haplo_proband_names
sample_annotation = c("GT", "AD", "DP", "GQ", "PL")
tidy_haplo_proband <- haplo_proband_list%>%
  rbindlist(use.names=TRUE, fill=TRUE, idcol="sample") %>%
  select(sample, Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, gnomAD_exome_ALL, SIFT_score, Otherinfo) %>%
  separate(Otherinfo, c("a", "b", "c", "d", "e", "f", "g", "h", "i", "filter_status", "k", "annotation", "annotation_results"), "\t") %>%
  select(-a, -b, -c, -d, -e, -f, -g, -h, -i, -k, -annotation) %>%
  separate(annotation_results, sample_annotation, ":")


# begin filtering tidy data frames ----------------------------------------


filter_criteria <- c('nonsynonymous SNV', 'stopgain', 'frameshift insertion', 'frameshift deletion')

laura_li_proband_counts <- tidy_laura_li_proband %>%
  filter(`Func.refGene` == 'exonic') %>%
  filter(ExonicFunc.refGene %in% filter_criteria) %>%
  filter(gnomAD_exome_ALL < 0.05) %>%
  group_by(`Chr`, `Start`) %>%
  mutate(variant_count = n())  %>%
  filter(variant_count > 0) %>%
  group_by(`Gene.refGene`, `ExonicFunc.refGene`) %>%
  mutate(gene_count = n_distinct(`Start`)) %>%
  arrange(desc(`gene_count`),`Chr`, `Start`)

# laura_li_parental_counts <- tidy_laura_li_parental %>%
#   filter(`Func.refGene` == 'exonic') %>%
#   filter(ExonicFunc.refGene %in% filter_criteria) %>%
#   filter(gnomAD_exome_ALL < 0.05) %>%
#   group_by(`Chr`, `Start`) %>%
#   mutate(variant_count = n())  %>%
#   filter(variant_count > 0) %>%
#   group_by(`Gene.refGene`, `ExonicFunc.refGene`) %>%
#   mutate(gene_count = n_distinct(`Start`)) %>%
#   arrange(desc(`gene_count`),`Chr`, `Start`)

haplo_proband_counts <- tidy_haplo_proband %>%
  filter(`Func.refGene` == 'exonic') %>%
  filter(ExonicFunc.refGene %in% filter_criteria) %>%
  filter(gnomAD_exome_ALL < 0.05) %>%
  group_by(`Chr`, `Start`) %>%
  mutate(variant_count = n())  %>%
  filter(variant_count > 0) %>%
  group_by(`Gene.refGene`, `ExonicFunc.refGene`) %>%
  mutate(gene_count = n_distinct(`Start`)) %>%
  arrange(desc(`gene_count`),`Chr`, `Start`)


# haplo_parental_counts <- tidy_haplo_parental %>%
#   filter(`Func.refGene` == 'exonic') %>%
#   filter(ExonicFunc.refGene %in% filter_criteria) %>%
#   filter(gnomAD_exome_ALL < 0.05) %>%
#   group_by(`Chr`, `Start`) %>%
#   mutate(variant_count = n())  %>%
#   filter(variant_count > 0) %>%
#   group_by(`Gene.refGene`, `ExonicFunc.refGene`) %>%
#   mutate(gene_count = n_distinct(`Start`)) %>%
#   arrange(desc(`gene_count`),`Chr`, `Start`)
  
sdf_ll_pro <- laura_li_proband_counts %>%
  select(`Gene.refGene`, `gene_count`, `gnomAD_exome_ALL`) %>%
  group_by(`Gene.refGene`) %>%
  summarise(num_variants_ll = n())
  

sdf_ll_par <- laura_li_parental_counts %>%
  group_by(`Gene.refGene`, `Chr`, `Start`, `End`, `gnomAD_exome_ALL`, `variant_count`) %>%
  summarise(n_avg_AF = mean(as.numeric(`AF`))) %>%
  filter(variant_count < 2)  

sdf_h_pro <- haplo_proband_counts %>%
  select(`Gene.refGene`, `gene_count`, `gnomAD_exome_ALL`) %>%
  group_by(`Gene.refGene`) %>%
  summarise(num_variants_haplo = n())
 

sdf_h_par <- haplo_parental_counts %>%
  group_by(`Gene.refGene`, `Chr`, `Start`, `End`, `gnomAD_exome_ALL`, `variant_count`) %>%
  summarise(n_avg_AF = mean(as.numeric(`AF`))) %>%
  filter(variant_count < 1)  

compdf <- sdf_ll_pro %>%
  inner_join(sdf_h_pro, by = c("Gene.refGene"))

compdf.plot = melt(compdf)

l <- density(compdf$num_variants_ll)
plot(l)
p <- ggplot(aes(x=value, color=variable), data=compdf.plot)
p + geom_density()

# haplotype comparisons ---------------------------------------------------





 

