#!/bin/bash

#This script is a pipeline for generating count tables from fastq files using kallisto.

usage()
{
    echo "usage: code to align FASTQ files to a kallisto index."
    echo "-i, --index: index"
    echo "-Q, --FASTQ: FASTQ file to grep"
    echo "-o, --output: output file, tab delimited count table"
    echo "Timothy J. Eisen, 2021 01 23"
}

main()
{

}

while [ "$1" != "" ]; do
    case $1 in
        -i | --index )          shift
                                index=$1
                                ;;
        -Q | --FASTQ )          shift
                                FASTQ=$1
                                ;;
        -o | --output )         shift
                                count_table=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done
main
