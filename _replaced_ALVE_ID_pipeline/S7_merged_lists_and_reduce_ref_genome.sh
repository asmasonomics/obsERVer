#!/bin/sh
## ./S7_merged_lists_and_reduce_ref_genome.sh prefix ref_genome.fa

# help message
if [ $# -ne 2 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./S7_merged_lists_and_reduce_ref_genome.sh prefix ref_genome.fa"; echo ""
   echo "Please rerun script with correct usage from the filtered bed file directory. The script will now exit."; echo ""
   exit
fi

## 1. Concatenate all filtered_softclipping.bed files and sort

cat *g.bed | sort -k1,1 -k2,2n > "${1}_all_positions_with_softclipping.bed"

## 2. Create a reduced reference genome file

awk '/^>/ {print (NR==1 ? "" :RS)$0; next}{printf "%s", $0}END{printf RS}' $2 > tmp_genome.fa
for i in `cut -f1 "${1}_all_positions_with_softclipping.bed" | sort | uniq`; do grep -A1 $i tmp_genome.fa >> "${1}_reduced_genome.fa"; done
rm tmp_genome.fa
samtools faidx "${1}_reduced_genome.fa"

