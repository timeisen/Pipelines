#!/bin/bash

#  RNAseq_pipeline_STAR_v6.sh
#
#   DOES NOT REMOVE SOFT CLIPPED READS
#  Created by Timothy Eisen on 6/19/15.
#  Updated by Timothy Eisen on 9/16/16.
#  Version 3 makes the following modifications:
#   changes FASTQ input scripts to accept .gz files
#   uses full path lengths
#   removes unnecessary intermediate files
#   creates a sorted bam file as the main output alignment
#   updates dependencies format
#   Version 4 adds proper prefix to the samtools sort function, to prevent overlapping files when running multiple jobs.
#   Version 6 is updated for the current version of the metabolic labeling paper. It outputs rpm, rpkm, and tags. 
#   Version 7 is an update for some neuro analyses. It's almost exactly the RNAseq pipeline except that it switches the intersect bed strandedness and doesn't remove random nucleotides.
#   Version 8 changes the first line to work with .gz, not .tar.gz, files. 

##Requires that at least one argument is present
if (! getopts ":1:s:b:d:" opt); then
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
while getopts ":1:s:b:i:d:h" opt; do
    case $opt in
        1)
            fastq1=$OPTARG
            ;;
        s)
            star_index=$OPTARG
            ;;
        b)
            bed_file=$OPTARG #sorted
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

##filename extraction
filename=$(basename ${fastq1%%.*})
extension="${fastq1: -7}" #Added this line and lines 94-98 on 2019 02 05 in order to allow .txt files to be fed into the pipeline. 
subex="${fastq1##*.}" #Added this line and lines 94-98 on 2019 02 05 in order to allow .txt files to be fed into the pipeline. 

##creates log file with input commands
echo "Log file $(date +%y_%m_%d)" > $directory/log_file.txt;
echo "These files were created with NEXTflex_pipeline_STAR_v7.sh" >>$directory/log_file.txt;
echo "Command line input: "$0" "$@>>$directory/log_file.txt;
echo "FASTQ1 file: $fastq1" >>$directory/log_file.txt;
echo "STAR index: $star_index" >>$directory/log_file.txt;
echo "BED file: $bed_file" >>$directory/log_file.txt;
echo "Output file directory: $directory" >>$directory/log_file.txt;

#$1 is the STAR index
#$2 is the fastq input, .gz file
#$3 is the bedfile
#$4 is the file directory in which to create intermediate files

#filename extraction
filename=$(basename ${fastq1%%.*})


#read filter
if [ "$subex" = "txt" ]; then
    bsub -q 18 -J "unzip_fastq1" "cat $fastq1>$directory/${filename}_fastq1.txt"
elif [ "$extension" = ".tar.gz" ]; then
    bsub -q 18 -J "unzip_fastq1" "tar xzfO $fastq1>$directory/${filename}_fastq1.txt"
elif [ "$subex" = "gz" ]; then
    bsub -q 18 -J "unzip_fastq1" "gunzip -c $fastq1>$directory/${filename}_fastq1.txt"
else
    echo 'cannot read fastq ext.'
    exit 1
fi

bsub -q 18 -J "sort" "sort -k1,1V -k2,2n $bed_file > $directory/${filename}_bed_sorted.txt" 

#bowtie commands, arguments
bsub -q 18 -w "ended(unzip_fastq1)" -J "STAR" "STAR --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --outReadsUnmapped Fastx --readFilesIn $directory/${filename}_fastq1.txt --outFileNamePrefix $directory/${filename}_ > $directory/${filename}_stdOut_logFile.txt"

#sam to bam conversion
bsub -q 18 -w "ended(STAR)" -J "sam_to_bam" "samtools view -bS $directory/${filename}_Aligned.out.sam | samtools sort -T $directory/${filename}_temp -@ 8 -O bam > $directory/${filename}.bam"

#bam to interesect bed output. This script takes a significant portion of the pipeline in terms of runtime. #S param flips the strand.
bsub -q 18 -w "ended(sam_to_bam) && ended(sort)" -J "bam_to_intersectBed" "intersectBed -bed -abam $directory/${filename}.bam -b $directory/${filename}_bed_sorted.txt -wb -S > $directory/${filename}_intersect_bed_output.txt"

#unique bed entries
bsub -q 18 -w "ended(bam_to_intersectBed)" -J "filter_intersect_bed" "python2.7 /lab/solexa_bartel/teisen/RNAseq/Scripts/general/intersect_bed_to_unique_entries.py $directory/${filename}_intersect_bed_output.txt $directory/${filename}_unique_intersect_bed.txt"

#computes expression
bsub -q 18 -w "ended(filter_intersect_bed)" -J "compute_expression" "cut -f 13,16,18 $directory/${filename}_unique_intersect_bed.txt | sort | uniq -c > $directory/${filename}_gene_assignments_compiled.txt"

#determines gene lengths based on bed file
bsub -q 18 -J "feature_length_determination" "python2.7 /lab/solexa_bartel/teisen/RNAseq/Scripts/compile_bed_lengths.py $bed_file $directory/bed_file_gene_lengths.txt"

#compiles both mapping data and determines rpkm
bsub -q 18 -w "ended(compute_expression)" -J "rpkm_compile" "Rscript /lab/solexa_bartel/teisen/RNAseq/Scripts/general/compute_expression_v1.R $directory/${filename}_gene_assignments_compiled.txt $directory/bed_file_gene_lengths.txt $directory/${filename}_Log.final.out $directory/${filename}_expression_values.txt $directory/${filename}_expression_values_10rpm_cutoff.txt"

#removes intermediate files
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files1" "rm -f $directory/${filename}_filtered.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files2" "rm -f $directory/${filename}_Aligned.out.sam"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files3" "rm -f $directory/${filename}_SJ.out.tab"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files4" "rm -f $directory/${filename}_Log.progress.out"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files5" "rm -f $directory/${filename}_bed_sorted.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files7" "rm -f $directory/bed_file_gene_lengths.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files8" "rm -f $directory/${filename}_stdOut_logFile.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files9" "rm -f $directory/${filename}_intersect_bed_output.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files10" "rm -f $directory/${filename}_unique_intersect_bed.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files11" "rm -f $directory/${filename}_gene_assignments_compiled.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files12" "rm -f $directory/${filename}_fastq1.txt"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files13" "rm -rf $directory/${filename}__STARtmp"

