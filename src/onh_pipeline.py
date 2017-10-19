#!/usr/bin/python

import subprocess
import sys
import re
import os
import datetime
import argparse
import gzip
import pipes
from ped_parser import Individual, Family

# list here all steps of the pipeline. Pipeline will run all these steps, if not requested otherwise in -s argument
steps_to_process_all = ["annovar", "vcf_anno"]


parser = argparse.ArgumentParser(description="runs onh pipeline")
parser.add_argument("-1", "--fastq-r1", dest="f1", help="fastq R1 file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-2", "--fastq-r2", dest="f2", help="fastq R2 file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-a", "--adapter", dest="adapter", help="comma separated list of adapter sequences to trim REQUIRED", metavar="ADAPTER", required=True)
parser.add_argument("-d", "--out", dest="out_dir", help="output directory REQUIRED", metavar="DIRECTORY", required=True)
parser.add_argument("-s", "--steps-to-run", dest="steps", help="steps of pipeline to run")
parser.add_argument("-rg", "--reference", dest="reference", default= "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", help="human reference genome (hg19/Grch38)", metavar="FILE[.fasta]", required=True)
parser.add_argument("--overwrite", dest="overwrite")

options = parser.parse_args()
fastq_r1_location = options.f1
fastq_r2_location = options.f2
reference_genome = options.reference
bwa_index = reference_genome.replace("WholeGenomeFasta/genome.fa", "BWAIndex/genome.fa")
output_directory = options.out_dir
adapter_sequence = options.adapter#.split(",")
output_directory_root = options.out_dir
if(options.steps == "All" or None):
	steps_to_process = steps_to_process_all
else:
	steps_to_process = options.steps.split(",")

print "Will run steps:", steps_to_process



num_threads = "7"
# *********************************************************************
# DEFINITION OF PATHS

Trimmomatic  = ["java", "-jar", "/usr/share/java/Trimmomatic-0.36/trimmomatic-0.36.jar"]
BwaMem = ["bwa", "mem"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java",  "-jar", "/usr/share/java/picard.jar", "SortSam"]
MarkDuplicates = ["java",  "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
Mosdepth = ["mosdepth"]
AddOrReplaceReadGroups = ["java",  "-jar", "/usr/share/java/picard.jar", "AddOrReplaceReadGroups"]
BuildBamIndex = ["java",  "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
BaseRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator"]
PrintReads = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "PrintReads"]
VariantFiltration = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "VariantFiltration"]
SelectVariants = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "SelectVariants"]
HaplotypeCaller = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
GenotypeGVCFs = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "GenotypeGVCFs"]
VariantRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "VariantRecalibrator"]
ApplyRecalibration = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "ApplyRecalibration"]
CalculateGenotypePosteriors = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "CalculateGenotypePosteriors"]
VariantAnnotator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "VariantAnnotator"]
TableAnnovar = ["/home/skevin/TOOLS/annovar/table_annovar.pl"]
VcfAnno = ["/home/skevin/TOOLS/bin/vcfanno_linux64"]
genmod_binary = ["genmod"]

min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
adapter_sequence = "AAGCAGTGGTATCAA"

r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
family_name = r1_base_filename[:3]

#~ family_vcf = gatk_output_dir+family_name+".vcf"
#~ if os.path.isfile(family_vcf):
	#~ sys.exit('this sample has already entered the pipeline')

# /DEFINITION OF PATHS
# *********************************************************************

# READ GROUP PARSING
# =====================================================================
#if("rg_parse" in steps_to_process) or ("all" in steps_to_process):
with gzip.open(fastq_r1_location) as fasta:
	s = fasta.readline().split(":")
	RGID =  "_".join(s[2:5])
	RGSM = r1_base_filename
	RGPL = "Illumina"
	RGLB = "lib"
	read_groups = '@RG\\tID:'+RGID+'\\tSM:'+RGSM+'\\tPL:'+RGPL
print read_groups

# TRIMMOMATIC
# =====================================================================
trimmomatic_output_dir = output_directory + "trimmomatic/"
fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + ".trimmed.fastq"
fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + ".trimmed.fastq"

if("trimmomatic" in steps_to_process):
	if not os.path(trimmomatic_output_dir):
		os.makedirs(trimmomatic_output_dir)
	print "running trimmomatic"
	cmd = trimmomatic_cmd[:]
	#print cmd
	cmd.append("PE")
	cmd.append("-threads")
	cmd.append(str(num_threads))
	#cmd.append("-trimlog")
	#cmd.append(log)
	cmd.append(fastq_r1_trimmed)
	cmd.append(fastq_r2_trimmed)
	cmd.append("LEADING:"+min_base_quality)
	cmd.append("TRAILING:"+min_base_quality)
	cmd.append("SLIDINGWINDOW:"+sliding_window)
	cmd.append("MINLEN:"+min_read_length)
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd


#~ # FASTQC
#~ # =====================================================================
#~ fastqc_output_dir = output_directory + "fastqc/"
#~ fastqc_output = fastqc_output_dir+r1_base_filename+"_fastqc.zip"
#~ if not os.path.isdir(fastqc_output_dir):
	#~ subprocess.call(["mkdir", fastqc_output_dir])
#~ if(not os.path.isfile(fastqc_output)) and(("fastqc" in steps_to_process) or ("all" in steps_to_process)):
	#~ cmd = ["fastqc",fastq_r1_location,fastq_r2_location]
	#~ print " ".join(cmd)
	#~ out_file = fastqc_output_dir + "fastqc.log"
	#~ err_file = fastqc_output_dir + "fastqc.err"
	#~ print " ".join(cmd)
	#~ with open(err_file,"w")as outerr:
		#~ with open(out_file,"w")as outf:
		 	#~ errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			#~ if(errcode == 0):
				#~ print "fastqc finished successfully"
			#~ else:
				#~ print "fastqc failed !!!!"
				#~ del steps_to_process[:]
	#~ del cmd

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 

bwa_sam = bwa_output_dir+r1_base_filename + ".sam"
print steps_to_process
if not os.path.isdir(bwa_output_dir):
	subprocess.call(["mkdir", bwa_output_dir])
if (not os.path.isfile(bwa_sam)) and(("bwa-aligner" in steps_to_process) or ("all" in steps_to_process)):
	print "beginning alignment"
	out_file = bwa_output_dir +r1_base_filename+".bwa.log"
	print fastq_r1, fastq_r2
	cmd = BwaMem[:]
	cmd.append("-M")
	cmd.append("-t")
	cmd.append("6")
	cmd.append("-R")
	cmd.append(read_groups)
	#cmd.append(log)
	cmd.append(bwa_index)
	cmd.append(fastq_r1_location) 
	cmd.append(fastq_r2_location)
	print " ".join(cmd)
        f = open(bwa_sam, "w")
        errcode = subprocess.call(cmd, stdout=f)
        if(errcode == 0):
                print "BWA finished successfully"
        else:
                print "BWA failed !!!!"
        f.close()
        del cmd



# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
picard_bam = picard_output_dir+r1_base_filename + ".bam"
if not os.path.isdir(picard_output_dir):
	subprocess.call(["mkdir", picard_output_dir])
if(not os.path.isfile(picard_bam)) and(("samformatconverter" in steps_to_process)):
	cmd = SamFormatConverter[:]
	cmd.append("I=")
	cmd.append(bwa_sam)
	cmd.append("O=")
	cmd.append(picard_bam)
	cmd.append("MAX_RECORDS_IN_RAM=200000")
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "samformatconverter finished successfully"
	else:
		print "samformatconverter failed !!!!"
		del steps_to_process[:]
	del cmd


# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_output_dir+r1_base_filename + "_sorted.bam"
if (not os.path.isfile(sortsam_bam)) and(("sortsam" in steps_to_process)):
	print "sortsam_bam:",sortsam_bam
	cmd = SortSam[:]
	cmd.append("I=")
	cmd.append(picard_bam)
	cmd.append("O=")
	cmd.append(sortsam_bam)
	cmd.append("SORT_ORDER=coordinate")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "sortsam finished successfully"
	else:
		print "sortsam failed !!!!"
		del steps_to_process[:]
	del cmd


# MARKDUPLICATES [and remove]
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
removed_duplicates_bam= picard_output_dir+ r1_base_filename + "_removed_duplicates.bam"
remove_dup_metrics = picard_output_dir+r1_base_filename+"_removed_dup_metrics.txt"
picard_tmp_dir = picard_output_dir+"tmp/"
if(not os.path.isfile(removed_duplicates_bam)) and (("removeduplicates" in steps_to_process)):
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O=")
	cmd.append(removed_duplicates_bam)
	cmd.append("M=")
	cmd.append(remove_dup_metrics)
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	cmd.append("REMOVE_DUPLICATES=true")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "markduplicates finished successfully"
	else:
		print "markduplicates failed !!!!"
		del steps_to_process[:]
	del cmd

# MOSDEPTH COVERAGE
# =====================================================================
#MOSDEPTH
coverage_table= picard_output_dir+r1_base_filename + "_coverage.txt"
coverage_dist= picard_output_dir+r1_base_filename + "_cumul.dist"
target_bed = "./bin/agilent_coverage_files/agilent_sureselect_v5_covered.sorted.bed"

if(not os.path.isfile(coverage_table) and "mosdepth" in steps_to_process):
	cmd = Mosdepth[:]
	cmd += ["--distribution", coverage_dist]
	cmd += ["--by", target_bed]
	cmd += [removed_duplicates_bam]
	print " ".join(cmd)
	log_file = picard_output_dir +r1_base_filename+"_mosdepth.log"
	with open(log_file,"w")as logf:
		with open(coverage_table, "w") as outf:
			errcode = subprocess.call(cmd, stdout = outf)
	del cmd

# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
remove_dup_index= picard_output_dir+r1_base_filename+ "_removed_duplicates.bam.bai"
if(not os.path.isfile(remove_dup_index)) and(("buildbamindex" in steps_to_process)):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(removed_duplicates_bam)
	cmd.append("O=")
	cmd.append(remove_dup_index)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "buildbamindex finished successfully"
	else:
		print "buildbamindex failed !!!!"
		del steps_to_process[:]
	del cmd

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
recal_report = gatk_output_dir+r1_base_filename+"_recal_report.table"
indels_vcf_dbsnp = "../Homo_sapiens/DBSNP/dbsnp_138.hg19.excluding_sites_after_129.vcf"
indels_vcf_Mills = "/dataVolume/storage/rb_pipeline/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "/dataVolume/storage/rb_pipeline/1000G_phase1.indels.hg19.sites.vcf"
if not os.path.isdir(gatk_output_dir):
	subprocess.call(["mkdir", gatk_output_dir])
if(not os.path.isfile(recal_report)) and (("baserecalibrator" in steps_to_process)):
	cmd = BaseRecalibrator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(removed_duplicates_bam)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_Mills)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_1000g)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-o")
	cmd.append(recal_report)
	cmd.append("-rf")
	cmd.append("BadCigar")
	out_file = gatk_output_dir + r1_base_filename+"_baserecalibrator.log"
	err_file = gatk_output_dir + r1_base_filename+"_baserecalibrator.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
			errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "baserecalibrator finished successfully"
			else:
				print "baserecalibrator failed !!!!"
	del cmd


# PRINTREADS
# =====================================================================
#PrintReads -T=PrintReads -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_bam = gatk_output_dir+r1_base_filename+"_recalibrated.bam"
if (not os.path.isfile(recal_bam)) and (("printreads" in steps_to_process)):
	cmd = PrintReads[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(removed_duplicates_bam)
	cmd.append("-BQSR")
	cmd.append(recal_report)
	cmd.append("-o")
	cmd.append(recal_bam)
	out_file = gatk_output_dir + r1_base_filename+"_printreads.log"
	err_file = gatk_output_dir + r1_base_filename+"_printreads.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "printreads finished successfully"
			else:
				print "printreads failed !!!!"
				del steps_to_process[:]
	del cmd
	
# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
recal_bam_index= gatk_output_dir+r1_base_filename+ "_recalibrated.bam.bai"
if(not os.path.isfile(recal_bam_index)) and (("recalbamindex" in steps_to_process)):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(recal_bam)
	cmd.append("O=")
	cmd.append(recal_bam_index)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "recalbamindex finished successfully"
	else:
		print "recalbamindex failed !!!!"
		del steps_to_process[:]
	del cmd

# G-HAPLOTYPECALLER
# =====================================================================
g_haplocaller_cell_dir = gatk_output_dir+family_name + "_gvcfs/"
haplocaller_gvcf = g_haplocaller_cell_dir+r1_base_filename+"_haplocaller.g.vcf"
if not os.path.isdir(g_haplocaller_cell_dir):
	subprocess.call(["mkdir", g_haplocaller_cell_dir])
if (not os.path.isfile(haplocaller_gvcf)) and ("g-haplotypecaller" in steps_to_process):
	cmd  = HaplotypeCaller[:]
	cmd += ["-R", reference_genome]
	cmd += ["-I", recal_bam]
	cmd += ["--dbsnp", indels_vcf_dbsnp]
	cmd += ["-L", target_intervals]
	cmd += ["-o", haplocaller_gvcf]
	cmd += ["-ERC", "GVCF"]
	cmd += ["-variant_index_type", "LINEAR"]
	cmd += ["-variant_index_parameter", "128000"]
	print ' '.join(str(p) for p in HaplotypeCaller)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "joint haplotypecaller for sample "+r1_base_filename+" finished successfully"
	else:
		print "joint haplotypecaller for sample "+r1_base_filename+" failed !!!!"
	del cmd

# GENOTYPEGVCFS
# =====================================================================
joint_vcf = gatk_output_dir+"onh_joint_call.vcf"
if(not os.path.isfile(joint_vcf)) and ("genotypegvcfs" in steps_to_process):
	cmd = GenotypeGVCFs[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	relpath = os.path.relpath(gatk_output_dir)
	for root, dirnames, filenames in os.walk(gatk_output_dir):
		for filename in filenames:
			if filename.endswith(".g.vcf"):
				filename = os.path.join(root, filename)
				cmd += ["--variant", filename]
	cmd.append("-stand_call_conf")
	cmd.append("50.0")
	cmd.append("-o")
	cmd.append(joint_vcf)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "genotypegvcfs for sample "+r1_base_filename+" finished successfully"
	else:
		print "genotypegvcfs for sample "+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd


# SNP VARIANTRECALIBRATOR
# =====================================================================
snps_omni = "./bin/1000G_omni2.5.hg19.sites.vcf"
snps_hapmap = "./bin/hapmap_3.3.hg19.sites.vcf"
snps_1000G = "./bin/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
snps_dbsnp = "../Homo_sapiens/DBSNP/dbsnp_138.hg19.excluding_sites_after_129.vcf"
vqsr_snp_recal_report = gatk_output_dir+"joint_recalibrate_SNP.recal"
vqsr_snp_tranches = gatk_output_dir+"joint_recalibrate_SNP.tranches"
vqsr_snp_plot_R = gatk_output_dir+"joint_recalibrate_SNP_plots.R"

if(not os.path.isfile(vqsr_snp_recal_report)) and ("snp_vqsr" in steps_to_process):
	cmd = VariantRecalibrator[:]
	cmd += ["-R", reference_genome]
	cmd += ["-input", joint_vcf]  
	cmd += ["-resource:hapmap,known=false,training=true,truth=true,prior=15.0", snps_hapmap]  
	cmd += ["-resource:omni,known=false,training=true,truth=true,prior=12.0", snps_omni]  
	cmd += ["-resource:1000G,known=false,training=true,truth=false,prior=10.0", snps_1000G]  
	cmd += ["-resource:dbsnp,known=true,training=false,truth=false,prior=2.0", snps_dbsnp]  
	cmd += ["-an", "DP"]
	cmd += ["-an", "QD"]
	cmd += ["-an", "FS"]
	cmd += ["-an", "SOR"]  
	cmd += ["-an", "MQ"] 
	cmd += ["-an", "MQRankSum"]  
	cmd += ["-an", "ReadPosRankSum"]  
	cmd += ["-mode", "SNP"]  
	cmd += ["-tranche", "100.0", "-tranche", "99.9", "-tranche", "99.0", "-tranche", "90.0"]  
	cmd += ["-recalFile", vqsr_snp_recal_report]  
	cmd += ["-tranchesFile", vqsr_snp_tranches]  
	cmd += ["-rscriptFile", vqsr_snp_plot_R] 
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "snp variant reclibrator for family "+family_name+" finished successfully"
	else:
		print "snp variant recalibrator for family "+family_name+" failed !!!!"
	del cmd

# SNP APPLY RECALIBRATION
# =====================================================================
denovo_vcf = gatk_output_dir+family_name+"_denovo.vcf"
snp_recalibrated_vcf = gatk_output_dir+"snp_recalibrated.vcf"
if(not os.path.isfile(snp_recalibrated_vcf)) and ("snp_apply_vqsr" in steps_to_process):
	cmd = ApplyRecalibration[:]
	cmd += ["-R", reference_genome ] 
	cmd += ["-input", joint_vcf]  
	cmd += ["-mode", "SNP"]
	cmd += ["--ts_filter_level", "99.0"  ]
	cmd += ["-recalFile", vqsr_snp_recal_report]
	cmd += ["-tranchesFile", vqsr_snp_tranches]  
	cmd += ["-o" ,snp_recalibrated_vcf ]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "snp apply reclibration for family "+family_name+" finished successfully"
	else:
		print "snp apply reclibrationannotation for family "+family_name+" failed !!!!"
	del cmd

# INDEL VARIANTRECALIBRATOR
# =====================================================================
indels_vcf_dbsnp = "../Homo_sapiens/DBSNP/dbsnp_138.hg19.excluding_sites_after_129.vcf"
indels_vcf_Mills = "/dataVolume/storage/rb_pipeline/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
vqsr_indel_recal_report = gatk_output_dir+"joint_recalibrate_indel.recal"
vqsr_indel_tranches = gatk_output_dir+"joint_recalibrate_indel.tranches"
vqsr_indel_plot_R = gatk_output_dir+"joint_recalibrate_indel_plots.R"

if(not os.path.isfile(vqsr_indel_recal_report)) and ("indel_vqsr" in steps_to_process):
	cmd = VariantRecalibrator[:]
	cmd += ["-R" ,reference_genome  ]
	cmd += ["-input" , snp_recalibrated_vcf ] 
	cmd += ["-resource:mills,known=false,training=true,truth=true,prior=12.0", indels_vcf_Mills]
	cmd += ["-resource:dbsnp,known=true,training=false,truth=false,prior=2.0", indels_vcf_dbsnp ]
	cmd += ["-an" ,"QD" ]
	cmd += ["-an", "DP" ] 
	cmd += ["-an", "FS"  ]
	cmd += ["-an" ,"SOR"  ]
	cmd += ["-an" ,"MQRankSum" ] 
	cmd += ["-an" ,"ReadPosRankSum"]  
	cmd += ["-an", "InbreedingCoeff"]
	cmd += ["-mode" ,"INDEL"  ]
	cmd += ["-tranche", "100.0", "-tranche", "99.9" ,"-tranche", "99.0", "-tranche", "90.0"]  
	cmd += ["--maxGaussians", "4"  ]
	cmd += ["-recalFile", "recalibrate_INDEL.recal" ] 
	cmd += ["-tranchesFile", "recalibrate_INDEL.tranches"]  
	cmd += ["-rscriptFile", "recalibrate_INDEL_plots.R" ]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "indel variant recalibrationannotation for family "+family_name+" finished successfully"
	else:
		print "indel variant recalibrationannotation for family "+family_name+" failed !!!!"
	del cmd
	
# INDEL APPLY RECALIBRATION
# =====================================================================
indel_recalibrated_vcf = gatk_output_dir+"indel_snp_recalibrated.vcf"
if(not os.path.isfile(indel_recalibrated_vcf)) and ("indel_apply_vqsr" in steps_to_process):
	cmd = ApplyRecalibration[:]
	cmd += ["-R" ,reference_genome]  
	cmd += ["-input", snp_recalibrated_vcf  ]
	cmd += ["-mode", "INDEL"  ]
	cmd += ["--ts_filter_level", "99.0"  ]
	cmd += ["-recalFile", "recalibrate_INDEL.recal"]  
	cmd += ["-tranchesFile", "recalibrate_INDEL.tranches" ] 
	cmd += ["-o", indel_recalibrated_vcf]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "indel apply recalibrationannotation for family "+family_name+" finished successfully"
	else:
		print "indel apply recalibrationannotation for family "+family_name+" failed !!!!"
	del cmd

# SELECT VARIANTS
# =====================================================================
family_vcf = gatk_output_dir+family_name+".vcf"
family_se = "^"+family_name+".*$"
if(not os.path.isfile(family_vcf)) and ("select_variants" in steps_to_process):
	cmd = SelectVariants[:]
	cmd += ["-R", reference_genome]
	cmd += ["-V", indel_recalibrated_vcf]
	cmd += ["-o", family_vcf]
	cmd += ["-se", family_se]
	cmd += ["-env", "-ef"]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "select_variants for sample "+r1_base_filename+" finished successfully"
	else:
		print "select_variants for sample "+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
    
# CALCULATE GENOTYPE POSTERIORS
# =====================================================================
cgp_vcf = gatk_output_dir+family_name+"_cgp.vcf"
joint_ped = gatk_output_dir+family_name+".ped"
if(not os.path.isfile(cgp_vcf)) and ("calc_gp" in steps_to_process):
	cmd = CalculateGenotypePosteriors[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd += ["--supporting", "./data/1000G_phase1.indels.hg19.sites.vcf"]
	cmd += ["-ped", joint_ped]
	cmd += ["-V", family_vcf]
	cmd += ["-o", cgp_vcf]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "calc_gp for family "+family_name+" finished successfully"
	else:
		print "calc_gp for family "+family_name+" failed !!!!"
		del steps_to_process[:]
	del cmd

# VARIANTFILTRATION
# =====================================================================
#VariantFiltration -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
filtered_vcf = gatk_output_dir+family_name+"_filtered.vcf"
if(not os.path.isfile(filtered_vcf)) and ("variantfiltration" in steps_to_process):
	cmd = VariantFiltration[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-V")
	cmd.append(cgp_vcf)
	cmd.append("--variant")
	cmd += ["-G_filter", "GQ < 20.0", "-G_filterName", "lowGQ"]
	cmd += ["-o", filtered_vcf]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "variantfiltration finished successfully"
	else:
		print "variantfiltration failed !!!!"
		del steps_to_process[:]
	del cmd

 # GATK ANNOTATE VARIANTS
# =====================================================================
annotated_vcf = gatk_output_dir+family_name+"_annotated.vcf"
if(not os.path.isfile(denovo_vcf)) and ("annotate_vars" in steps_to_process):
	cmd = VariantAnnotator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd += ["-V", filtered_vcf]
	cmd += ["-A", "PossibleDeNovo"]
	cmd += ["-A", "VariantType"]
	cmd += ["-ped", joint_ped]
	cmd += ["--dbsnp", indels_vcf_dbsnp]
	cmd += ["-o", annotated_vcf]
	
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "de novo annotation for family "+family_name+" finished successfully"
	else:
		print "de novo annotation for family "+family_name+" failed !!!!"
	del cmd

# ANNOVAR
# =====================================================================
annovar_output_dir = output_directory+"annovar/"
annovar_vcf = annovar_output_dir+family_name+"_annovar"

if not os.path.isdir(annovar_output_dir):
	subprocess.call(["mkdir", annovar_output_dir])
if(not os.path.isfile(annovar_vcf+".hg19_multianno.vcf")) and (("annovar" in steps_to_process)):

# table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene -operation g -nastring .
	cmd = TableAnnovar[:]
	cmd.append(annotated_vcf)
	cmd.append("/home/skevin/TOOLS/annovar/humandb/")
	cmd.append("-buildver")
	cmd.append("hg19")
	cmd.append("--out")
	cmd.append(annovar_vcf)
	cmd.append("-remove")
	cmd.append("-protocol")
	cmd += ["refGene"]
	cmd.append("-operation")
	cmd += ["g"]
	cmd.append("-nastring")
	cmd.append(".")
	cmd.append("-otherinfo")
	cmd.append("-vcfinput")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "annovar finished successfully"
	else:
		print "annovar failed !!!!"
		del steps_to_process[:]
	del cmd

# VCF ANNO
# =====================================================================
anno_vcf = gatk_output_dir+family_name+"_anno.vcf"
if(not os.path.isfile(anno_vcf)) and ("vcf_anno" in steps_to_process):
	cmd = VcfAnno[:]
	cmd += ["-lua", "./bin/gnomad.lua"]
	cmd += ["./bin/gnomad.toml"]
	cmd += [annovar_vcf+".hg19_multianno.vcf"]
	print " ".join(cmd)
	log_file = gatk_output_dir +family_name+"_vcf_anno.log"
	with open(log_file,"w")as logf:
		with open(anno_vcf, "w") as outf:
			errcode = subprocess.call(cmd, stdout = outf)
		del cmd

#~ # GENMOD
#~ # =====================================================================
genmod_vcf = gatk_output_dir+family_name+"_genmod.vcf"
ped_file = gatk_output_dir + family_name + ".ped"
if(not os.path.isfile(genmod_vcf)) and ("genmod" in steps_to_process):
	cmd = genmod_binary[:]
	cmd.append("annotate")
	cmd.append(anno_vcf)
	cmd += ["--annotate_regions"]
	print " ".join(cmd)
	print datetime.datetime.today()
	print ped_file
	log_file = gatk_output_dir +family_name+"_genmod.log"
	with open(log_file,"w")as logf:
		with open("genmod.temp", "w") as outf:
			annotate = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=logf)
		 	models = subprocess.Popen(["genmod", "models", "-", "--family_file", ped_file, "-t", "ped"], stdin=annotate.stdout, stdout = outf, stderr=logf)
		 	annotate.stdout.close()
			annotate.wait()
		 	models.communicate()
	del cmd
	
	with open(genmod_vcf, "w") as outf:
		cmd = ["awk", '$8~/GeneticModels/ || $1~/#/', "genmod.temp"]
		print " ".join(cmd)
		errcode = subprocess.call(cmd, stdout = outf)
		
		print datetime.datetime.today()
		del cmd

