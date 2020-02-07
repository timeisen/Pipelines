#!/bin/bash

#  RNAseq_pipeline_STAR_v8.sh
#
#	DOES NOT REMOVE SOFT CLIPPED READS
#  Created by Timothy Eisen on 6/19/15.
#  Updated by Timothy Eisen on 9/16/16.
#  Version 3 makes the following modifications:
#	changes FASTQ input scripts to accept .gz files
#	uses full path lengths
#	removes unnecessary intermediate files
#	creates a sorted bam file as the main output alignment
#   updates dependencies format
#   Version 4 adds proper prefix to the samtools sort function, to prevent overlapping files when running multiple jobs.
#   Version 6 is updated for the current version of the metabolic labeling paper. It outputs rpm, rpkm, and tags. 
#	Version 7 is an update for some neuro analyses.
#   Version 8 changes the fastq file parser to accept .gz rather than .tar.gz input fastq files. 
#   Version 9 uses the featureCounts software 2020 01 16, no standards, doesn't use a bed file, uses cutadapt. 

##Requires that at least one argument is present
if (! getopts ":1:s:g:d:" opt); then
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` \nRequired arguments: \n-1 Fastq1 \n-s STAR index \n-b Sorted bed annotation file \n-d Directory\n";
    exit $E_OPTERROR;
fi

##requires that the correct number of arguments is present
if [ "$#" -ne 8 ]; then
    echo -e "\nIllegal number of arguments (Error #2)"
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` \nRequired arguments: \n-1 Fastq1 \n-s STAR index \n-b Sorted bed annotation file\n-d Directory\n";
    exit
fi

##parses argument flags
while getopts ":1:s:g:i:d:h" opt; do
    case $opt in
        1)
            fastq1=$OPTARG
            ;;
        s)
            star_index=$OPTARG
            ;;
        g)
            gtf_file=$OPTARG #sorted
            ;;
   
        d)
            output_directory=$OPTARG
            ;;
        h)
            echo -e "-h print this screen and exit"
            exit $E_OPTERROR
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

directory=$output_directory
echo "$directory used as the working directory"

##creates directory in which to store intermediate and output files
mkdir $directory

##creates log file with input commands
echo "Log file $(date +%y_%m_%d)" > $directory/log_file.txt;
echo "These files were created with RNAseq_pipeline_STAR_v9.sh" >>$directory/log_file.txt;
echo "Command line input: "$0" "$@>>$directory/log_file.txt;
echo "FASTQ1 file: $fastq1" >>$directory/log_file.txt;
echo "STAR index: $star_index" >>$directory/log_file.txt;
echo "GTF file: $gtf_file" >>$directory/log_file.txt;
echo "Output file directory: $directory" >>$directory/log_file.txt;

#$1 is the STAR index
#$2 is the fastq input, .gz file
#$3 is the gtf_file
#$4 is the file directory in which to create intermediate files

#filename extraction
filename=$(basename ${fastq1%%.*})
extension="${fastq1: -7}" #Added this line and lines 94-98 on 2019 02 05 in order to allow .txt files to be fed into the pipeline. 
subex="${fastq1##*.}" #Added this line and lines 94-98 on 2019 02 05 in order to allow .txt files to be fed into the pipeline. 

#read filter
if [ "$subex" = "txt" ]; then
    bsub -q 18 -J "unzip_fastq1" "cat $fastq1 > $directory/${filename}_fastq1.txt"
elif [ "$extension" = ".tar.gz" ]; then
    bsub -q 18 -J "unzip_fastq1" "tar xzfO $fastq1 > $directory/${filename}_fastq1.txt"
elif [ "$subex" = "gz" ]; then
    bsub -q 18 -J "unzip_fastq1" "gunzip -c $fastq1 > $directory/${filename}_fastq1.txt"
else
    echo 'cannot read fastq ext.'
    exit 1
fi

#read filter, discards reads < 10 nt.
bsub -q 18 -w "ended(unzip_fastq1)" -J "cutadapt" "cutadapt -j 10 $directory/${filename}_fastq1.txt -u 8 -m 10 -a TCGTATGCCGTCTTCTGCTTG > $directory/${filename}_trimmed_fastq.txt"

#bowtie commands, arguments
bsub  -q 18 -w "ended(cutadapt)" -J "STAR" "STAR --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --readFilesIn $directory/${filename}_trimmed_fastq.txt --outFileNamePrefix $directory/${filename}_ > $directory/${filename}_stdOut_logFile.txt"

#index the bam file
bsub  -q 18 -w "ended(STAR)" -J "index" "samtools index $directory/${filename}_Aligned.sortedByCoord.out.bam"

#computes expression
bsub  -q 18 -w "ended(index)" -J "compute_expression" "featureCounts -T 10 -s 1 -f -a $gtf_file -o $directory/${filename}_featureCounts.txt $directory/${filename}_Aligned.sortedByCoord.out.bam"

#Removes raw fastq files. 
bsub -q 18 -w "ended(STAR)" -J "cleanup" "rm $directory/${filename}_fastq1.txt $directory/${filename}_trimmed_fastq.txt"

