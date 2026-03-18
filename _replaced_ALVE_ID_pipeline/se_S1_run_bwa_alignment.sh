#!/bin/sh
## run_bwa_alignment.sh path/to/read.fastq run_number lane_number genome/to/map.fa &> line_name.log

# help message
if [ $# -ne 4 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./run_bwa_alignment.sh path/to/read.fastq run_number lane_number genome/to/map.fa"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

outname=`echo $1 | rev | cut -d'/' -f1 | rev | sed s/.fastq.gz//g`

bwa mem -t 4 -R "@RG\tID:${outname}_$2_$3\tLB:${outname}\tPL:ILLUMINA\tSM:${outname}\tPU:$2_$3" $4 $1 > ${outname}.sam
samtools view -bhS ${outname}.sam > ${outname}.bam
samtools sort -@ 4 -m 4G ${outname}.bam ${outname}.sorted
samtools index ${outname}.sorted.bam

rm ${outname}.sam ${outname}.bam
