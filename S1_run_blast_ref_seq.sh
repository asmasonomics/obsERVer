# run_blast_ref_seq.sh ref_seq.fa ref_genome_blast_db_root

#!/bin/sh

# help message
if [ $# -ne 2 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./run_blast_ref_seq.sh ref_seq.fa ref_genome_blast_db_root"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

# split reference file into individual files
python ./accessory_scripts/seq_formatter.py $1 &> /dev/null

# perform a BLASTn search for each reference sequence and modify results file
for ref in `ls | grep ".fas$"`
do
    blastn -db $2 -query $ref -evalue 0.0000000001 -num_threads 8 -out $ref.results -outfmt "6 sacc sstart send" -max_target_seqs 1000000
    refname=`echo $ref | cut -d'_' -f1 | sed 's/-/_/g'`
    awk -v var="$refname" '{if ($2>$3){print $1 "\t" $3 "\t" $2 "\t-\t" var}else{print $1 "\t" $2 "\t" $3 "\t+\t" var}}' $ref.results | sort -k1,1 -k2,2n > $ref.out
    rm $ref.results 
    rm $ref
done

# concatenate all results files
cat *.fas.out | sort -k1,1 -k2,2n | uniq > hom_pos.txt
rm *.fas.out

# merge positions
python ./accessory_scripts/pos_merger.py hom_pos.txt &> /dev/null
rm hom_pos.txt
mv merged_positions.txt homology_positions.txt

# create positions files for each type of element identified, grouped by clusters (due to high homology and often multiple hits) 
grep -E "EAV|evJ|ART" homology_positions.txt | grep -E "EAV_HP|evJ" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tEAV-HP"}' > EAV-HP_positions.txt
grep -E "EAV|evJ|ART" homology_positions.txt | grep -v -E "EAV_HP|evJ" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tEAV"}'  > EAV-ART_positions.txt
grep "IC10" homology_positions.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tIC10"}'  > IC10_positions.txt
grep "rav0_5prime" homology_positions.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\trav0-5prime"}' > rav0-5prime_positions.txt
grep -E "RSV|ALV_A|ALV_L|ALVE|ev21|ev_C|ev3|ev6|ev1" homology_positions.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tALVE"}' > ALVE_positions.txt

cat EAV-HP_positions.txt EAV-ART_positions.txt IC10_positions.txt rav0-5prime_positions.txt ALVE_positions.txt | sort -k1,1 -k2,2n > homology_positions_with_matching_seq.txt
mv homology_positions.txt homology_positions_with_original_matches.txt

# perform sequence extraction
python ./accessory_scripts/seq_extract.py --prefix EAV-HP EAV-HP_positions.txt $2.fa &> /dev/null
python seq_extract.py --prefix EAV-ART EAV-ART_positions.txt $2.fa &> /dev/null
python ./accessory_scripts/seq_extract.py --prefix IC10 IC10_positions.txt $2.fa &> /dev/null
python ./accessory_scripts/seq_extract.py --prefix rav0-5prime rav0-5prime_positions.txt $2.fa &> /dev/null
python ./accessory_scripts/seq_extract.py --prefix ALVE ALVE_positions.txt $2.fa &> /dev/null
rm *positions.txt
