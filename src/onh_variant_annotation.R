#!/usr/bin/Rsript


# load required libraries -------------------------------------------------


library(VariantAnnotation)
library(biobroom)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tibble)
library(dplyr)
library(data.table)
library(purrr)

# load input --------------------------------------------------------------

# lli_filenames <- list.files(path="./laurali_output/ll_gno_anno_vcfs", pattern="*.vcf", full.names = TRUE)
annotated_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]_anno.vcf.gz$", full.names = TRUE)
genmod_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]_genmod.vcf.gz$", full.names = TRUE)

# lli_names <- substr(lli_filenames,35,41)
annotated_names <-substr(annotated_filenames,15,17) #format for haplocaller samples
genmod_names <-substr(genmod_filenames,15,17) #format for haplocaller samples

# read in vcfs ------------------------------------------------------------

# lli_list <- map(lli_filenames, readVcf, "hg19")
annotated_list <- lapply(annotated_filenames, function(x)readVcf(x, "hg19"))
genmod_list <- lapply(genmod_filenames, function(x)readVcf(x, "hg19"))

# names(lli_list) <- lli_names
names(annotated_list) <- annotated_names
names(genmod_list) <- genmod_names

anno_vcf <-  annotated_list[[21]]
genmod_vcf <- genmod_list[[21]]
# lli_vcf <- lli_list[[1]]
# lli_evcf <- S4Vectors::expand(lli_vcf)

#genmod_list <-  genmod_list[1:2]


# define required functions -----------------------------------------------

 # single_lli <- function(lli_vcf){
 # 
 #    lli_evcf <- S4Vectors::expand(lli_vcf)  
 # 
 #    anno_df <- data.frame(rowRanges(lli_evcf),
 #                          GT=geno(lli_evcf)$GT,
 #                          DP=info(lli_evcf)$DP,
 #                          AF=info(lli_evcf)$AF,
 #                          FC=as.character(info(lli_evcf)$FC),
 #                          GENEID=as.character(info(lli_evcf)$SGGI),
 #                          TXID=as.character(info(lli_evcf)$SGTI),
 #                          gnomad.AF=info(lli_evcf)$gno_af_all,
 #                          gno_filter=info(lli_evcf)$gno_filter,
 #                          gno_id=info(lli_evcf)$gno_id
 #    ) %>% 
 #      dplyr::rename_all(~gsub('^.*.pjt$', 'GT', .))
 #    
 #    ###################################################
 #    ### code chunk number 30: predictCoding_frameshift
 #    ###################################################
 #    ## CONSEQUENCE is 'frameshift' where translation is not possible
 #    
 #    
 #    anno_df <- dplyr::filter(anno_df, !is.na(GENEID))
 #    
 #  }
  
single_tidy <- function(anno_vcf, genmod_vcf){
  
  anno_evcf <- S4Vectors::expand(anno_vcf)  
  genmod_evcf <- S4Vectors::expand(genmod_vcf)
  
  
  anno_df <- data.frame(rowRanges(anno_evcf),
                        snp_id=info(anno_evcf)$dbsnp_id,
                        GT=geno(anno_evcf)$GT,
                        AF=info(anno_evcf)$AF, 
                        QD=info(anno_evcf)$QD,
                        AD=geno(anno_evcf)$AD,
                        DP=geno(anno_evcf)$DP,
                        FS=info(anno_evcf)$FS,
                        MQ=info(anno_evcf)$MQ,
                        CAF=unstrsplit(info(anno_evcf)[,"dbsnp_af"], " "),
                        VQSLOD=info(anno_evcf)$VQSLOD,
                        gnomad.wes.AF=info(anno_evcf)$gno_wes_af_all,
                        gnomad.wgs.AF=info(anno_evcf)$gno_wgs_af_all,
                        gno_filter=info(anno_evcf)$gno_filter,
                        gno_id=info(anno_evcf)$gno_id,
                        SOR=info(anno_evcf)$SOR,
                        MQRankSum=info(anno_evcf)$MQRankSum,
                        ReadPosRankSum=info(anno_evcf)$ReadPosRankSum,
                        VT=info(anno_evcf)$VariantType,
                        Func.refGene=info(anno_evcf)$Func.refGene,
                        Gene.refGene=unstrsplit(info(anno_evcf)$Gene.refGene),
                        GeneDetail.refGene=unstrsplit(info(anno_evcf)$GeneDetail.refGene),
                        ExonicFunc.refGene=unstrsplit(info(anno_evcf)$ExonicFunc.refGene),
                        AAChange.refGene=unstrsplit(info(anno_evcf)$AAChange.refGene),
                        hiDN=info(anno_evcf)$hiConfDeNovo,
                        loDN=info(anno_evcf)$loConfDeNovo
  ) %>% 
    dplyr::rename_all(~gsub('\\.\\d+\\.', '\\.', .)) %>% 
    dplyr::rename_all(~gsub('\\.\\d+_', '_pro_', .)) 
  
  variants <- group_by(anno_df, seqnames, start, end) %>% 
    filter(row_number() == 1) %>% 
    ungroup()

  variants$CAF <- sapply(strsplit(as.character(variants$CAF), " "), "[", 2)
  variants$CAF <- as.numeric(variants$CAF)
  ###################################################
  ### code chunk number 32: genetic models
  ###################################################  
  
  genmod_df <- data.frame(rowRanges(genmod_evcf),
                          GeneticModels=unlist(info(genmod_evcf)$GeneticModels)) %>%  
    dplyr::select(seqnames, start, end, REF, ALT, GeneticModels) %>% 
    mutate(GeneticModels = gsub("^.*:", "", GeneticModels)) 
  
  join_df <- left_join(variants, genmod_df, by=c("seqnames", "start", "end", "REF", "ALT")) %>% 
    dplyr::filter(!is.na(Gene.refGene)) %>% 
    dplyr::filter(ExonicFunc.refGene != "synonymous_SNV") %>% 
    dplyr::group_by(snp_id) %>% 
    dplyr::filter(row_number() == 1)
    
}

# test <- single_tidy(anno_vcf, genmod_vcf)

collate_vcfs <- function(anno_vcf_list, genmod_vcf_list){


  evcf_list <- mapply(single_tidy, annotated_list, genmod_list, SIMPLIFY = FALSE)
  tidy_vcfs <- data.table::rbindlist(evcf_list, idcol = "sample", fill = TRUE)
  
 library(rfPred)
  
  rfp_input <- dplyr::select(data.frame(tidy_vcfs), chr = seqnames, pos = start, ref = REF, alt = ALT) %>% 
    mutate(chr = gsub("chr", "", chr))
  
  rfp0 <- rfPred_scores(variant_list=rfp_input,
                        data="./bin/all_chr_rfPred.txtz",
                        index="./bin/all_chr_rfPred.txtz.tbi")
  
  tidy_vcfs <- mutate(tidy_vcfs, seqnames = gsub("chr","", seqnames))
  
  rfp <- left_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration")) %>% 
    mutate(VAF.Dad = AD.Dad_1.2/DP.Dad_1) %>% 
    mutate(VAF.Mom = AD.Mom_1.2/DP.Mom_1) %>% 
    mutate(VAF.pro = AD_pro_1.2/DP_pro_1) 
}


tidy_vcfs0 <- collate_vcfs(annotated_list, genmod_list)
saveRDS(tidy_vcfs0, "./results/tidy_vcfs0_20171004.rds")

tidy_vcfs0 <- readRDS("./results/tidy_vcfs0_20171004.rds")

retidy_vcfs <- function(my_vcfs){
  browser()
  variants <- dplyr::filter(my_vcfs, !is.na(Gene.refGene)) %>% 
    dplyr::filter(ExonicFunc.refGene != "synonymous_SNV") %>% 
    dplyr::group_by(sample, snp_id) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::group_by(snp_id) %>% 
    dplyr::filter(AD_pro_1.1 != 0 | AD_pro_1.2 != 0) %>% 
    dplyr::filter(AD.Dad_1.1 != 0 | AD.Dad_1.2 != 0) %>% 
    dplyr::filter(AD.Mom_1.1 != 0 | AD.Mom_1.2 != 0) %>%
    filter(GT.Dad_1 != "0/0" & GT.Mom_1 != "0/0" & GT_pro_1 != "0/0") %>% 
    dplyr::filter(gnomad.wes.AF < 0.10 | is.na(gnomad.wes.AF)) %>% 
    dplyr::filter(gnomad.wgs.AF < 0.10 | is.na(gnomad.wgs.AF)) %>% 
    dplyr::filter(CAF < 0.10 | is.na(CAF)) %>% 
    dplyr::filter(!is.na(sample)) %>% # check up on significance of this threshold
    group_by(snp_id) %>% 
    mutate(recurrence=paste(sample,collapse=';')) %>% 
    mutate(counts = n()) %>%    
    dplyr::arrange(desc(counts)) 
  
  genes <- variants %>%
    group_by(snp_id) %>% 
    filter(row_number() == 1) %>% 
    group_by(Gene.refGene) %>%
    mutate(gene_counts = sum(counts)) %>%
    mutate(gene_recurrence = paste(recurrence, collapse=";")) %>% 
    mutate(gene_recurrence = map_chr(strsplit(as.character(gene_recurrence) ,";"), function(x) paste(unique(x), collapse=";"))) %>% 
    mutate(gene_recurrence_counts = map_int(strsplit(as.character(gene_recurrence) ,";"), function(x) length(unique(x)))) %>%
    dplyr::arrange(desc(counts)) %>%
    ungroup() 
  
  samples <- variants %>% 
    dplyr::group_by(sample, snp_id) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>% 
    dplyr::select(sample, Gene.refGene) %>%
    group_by(sample) %>% 
    dplyr::summarise(Gene.refGene = paste(Gene.refGene, collapse=";")) %>% 
    dplyr::arrange(desc(sample)) 
  
  
  my_list <-  list("variants" = variants, "genes" = genes, "samples" = samples)
  
}


tidy_vcfs <- retidy_vcfs(tidy_vcfs0) 
# dplyr::select(-GT.394MOM_1, -AD.394MOM_1.1, -AD.394MOM_1.2, -DP.394MOM_1) 

tidy_base_path = "./results/20171004_tidy_vcfs/"
dir.create(tidy_base_path)
variants_path <- paste(tidy_base_path, "variants", "_tidy_table.csv", sep = "")
genes_path <- paste(tidy_base_path, "genes", "_tidy_table.csv", sep = "")
samples_path <- paste(tidy_base_path, "samples", "_tidy_table.csv", sep = "")

write.table(tidy_vcfs$variants, variants_path, sep = ",", row.names = FALSE)
write.table(tidy_vcfs$genes, genes_path, sep = ",", row.names = FALSE)
write.table(tidy_vcfs$samples, samples_path, sep = ",", row.names = FALSE)

my_path <- (c(variants_path, genes_path, samples_path))
tidy_vcfs <- lapply(my_path, read.table, sep=",", header  = TRUE)
tidy_vcfs <- setNames(tidy_vcfs, c("variants", "genes", "samples"))

