#!/bin/bash

# argv[1] = file with fastq to process
# argv[2] = organism
# argv[3] = output dir
# argv[4] = step to process

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
 	argv[1] = file with fastq to process
	argv[2] = organism
	argv[3] = output dir
	argv[4] = step to process'
    exit 0
fi

list_to_process=$1 #"input_fastq"
organism=$2
outdir=$3
step_to_process=$4

for b in `cat $list_to_process`;do 
	j=`echo $b | sed 's/_1/_2/'`
	cell_name=`echo $b | sed 's#.*/\(.*\)_1.*#\1#'`
	echo $b $j $cell_name
	#./single_cell_pipeline.py $i $j $organism $outdir $step_to_process
	./src/onh_pipeline.py -1 $b -2 $j \
	-a AAGCAGTGGTATCAA -d $outdir \
	-rg "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" \
	-s $step_to_process 
done
