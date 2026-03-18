#!/bin/sh
## extract_ref_seq_mapped_reads.sh file.bam

# help message
if [ $# -ne 1 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./extract_ref_seq_mapped_reads.sh file.bam"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

outname=`echo $1 | cut -d. -f1`

## 1. Extract header names of reads in given BAM file which map to pseudochromosome

samtools view $1 | awk '$3=="pseudochromosome"{print $1}' > "${outname}_matched_headers.tmp"
sort "${outname}_matched_headers.tmp" | uniq > "${outname}_matched_headers.txt"
rm "${outname}_matched_headers.tmp"

## 2. Extract BAM file headers and all matched reads from part 1 with their pairs, and sort the output

samtools view -H $1 > "${outname}_headers.sam"
samtools view $1 | grep -F -f "${outname}_matched_headers.txt" | awk '$2<200{print $0}' > "${outname}_reads.sam"
cat "${outname}_headers.sam" "${outname}_reads.sam" > "${outname}_viral_mapped.sam" 
samtools view -bhS "${outname}_viral_mapped.sam" > "${outname}_viral_mapped.bam"
samtools sort -n "${outname}_viral_mapped.bam" "${outname}_viral_mapped_sorted"
rm "${outname}_matched_headers.txt" "${outname}_headers.sam" "${outname}_reads.sam" "${outname}_viral_mapped.bam" "${outname}_viral_mapped.sam"

## 3. Extract and zip reads into fastq format 

samtools view "${outname}_viral_mapped_sorted.bam" | awk '{print "@"$1"#0/1\n"$10"\n+\n"$11}' > "${outname}_viral_mapped.fastq"; gzip "${outname}_viral_mapped.fastq"
rm "${outname}_viral_mapped_sorted.bam"


