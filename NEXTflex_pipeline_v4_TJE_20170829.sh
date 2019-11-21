#!/bin/bash

#  TAILseq_pipeline_v1_TJE.sh
#
#
#  Created by Timothy Eisen on 6/19/15.
#  Updated by Timothy Eisen on 10/03/15.
#  Updated by Timothy Eisen on 1/30/16 to create TAIL-seq modifications
#  Updated by Timothy Eisen on 2/22/16
#  Updated by Timothy Eisen on 6/1/17
#  Current version makes the following changes:
#  Does not remove soft clipped reads
#  #V4 calculates RPM, RPKM automatically

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

directory=$output_directory""_$(date +%y_%m_%d)
echo "$directory used as the working directory"

##creates directory in which to store intermediate and output files
mkdir $directory

##filename extraction
filename=$(basename ${fastq1%%.*})

##creates log file with input commands
echo "Log file $(date +%y_%m_%d)" > $directory/log_file.txt;
echo "Command line input: "$0" "$@>>$directory/log_file.txt;
echo "FASTQ1 file: $fastq1" >>$directory/log_file.txt;
echo "STAR index: $star_index" >>$directory/log_file.txt;
echo "BED file: $bed_file" >>$directory/log_file.txt;
echo "Output file directory: $directory" >>$directory/log_file.txt;

##unzip data
bsub -q 18 -J "unzip_fastq1" "tar xzfO $fastq1>$directory/${filename}_fastq1.txt"

##rev comp the fastq file
bsub -q 18 -w "ended(unzip_fastq1)" -J "revcomp" "fastx_reverse_complement -i $directory/${filename}_fastq1.txt -o $directory/${filename}_fastq_revcomp.txt"

#uncomment this line to use unzipped data, also comment out the above lines. 
#bsub -q 18 -J "revcomp" "fastx_reverse_complement -i $fastq1 -o $directory/${filename}_fastq_revcomp.txt"

##STAR commands, arguments
bsub -q 18 -w "ended(revcomp)" -J "STAR" "STAR --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --outReadsUnmapped Fastx --readFilesIn $directory/${filename}_fastq_revcomp.txt --outFileNamePrefix $directory/${filename}_map1_ > $directory/${filename}_stdOut_logFile.txt"

##sam to bam conversion
bsub -q 18 -w "ended(STAR)" -J "sam_to_bam" "samtools view -bS $directory/${filename}_map1_Aligned.out.sam | samtools sort -T $directory/${filename}_temp -@ 8 -O bam > $directory/${filename}_map1.bam"

##bam to interesect bed output. This script takes a significant portion of the pipeline in terms of runtime.
bsub -q 18 -w "ended(sam_to_bam)" -J "bam_to_intersectBed" "intersectBed -bed -abam $directory/${filename}_map1.bam -b $bed_file -wo -s > $directory/${filename}_wo_output.txt"

bsub -q 18 -w "ended(bam_to_intersectBed)" -J "compile" "cut -f 16 $directory/${filename}_wo_output.txt | sort | uniq -c | sort -rn > $directory/${filename}_gene_assignment_compiled.txt"

#compiles both mapping data and determines rpkm
bsub -q 18 -J "feature_length_determination" "python /lab/solexa_bartel/teisen/RNAseq/Scripts/compile_bed_lengths.py $bed_file $directory/bed_file_gene_lengths.txt"

bsub -q 18 -w "ended(compile) && ended(feature_length_determination)" -J "rpkm_compile" "Rscript /lab/solexa_bartel/teisen/RNAseq/Scripts/general/compute_expression_v2.R $directory/${filename}_gene_assignment_compiled.txt $directory/bed_file_gene_lengths.txt $directory/${filename}_map1_Log.final.out $directory/${filename}_expression_values.txt $directory/${filename}_expression_values_10rpm_cutoff.txt"

#removes intermediate files
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files1" "find $directory/ !  \( -name '${filename}_gene_assignment_compiled.txt' -o -name 'log_file.txt' -o -name '${filename}_map1_Log.final.out' -o -name '${filename}_expression_values.txt' -o -name '${filename}_expression_values_10rpm_cutoff.txt' \) -type f -exec rm -r {} +"
# bsub -q 18 -w "ended(rpkm_compile)" -J "remove_intermediate_files2" "rm -rf $directory/${filename}_map1__STARtmp/"



