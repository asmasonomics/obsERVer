## python extract_putative_sites.py [-h] viral_mapped.sorted.bam ref_homology_pos.bed ref_genome.fa

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(
	description="extract_putative_sites.py takes the viral_mapped bam file, converts to bed and then extracts putative insert regions",
	epilog="Author: Andrew Mason; Release: 19/04/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("bam_file", help="specify viral mapped BAM file for analysis")
parser.add_argument("ref_pos", help="reference positions identified by homology search (BED)")
#parser.add_argument("ref_gen", help="reference genome with sequences on individual lines (run fasta1line if necessary)")
usr_args = parser.parse_args()

file_prefix = ((usr_args.bam_file).split("/")[-1]).split(".")[0]
line_prefix = (((usr_args.bam_file).split("_viral")[0]).split("_")[-1]).replace("line","").replace("indiv","")
for res_file in [(file_prefix + "_filtered_locations.sam"),(file_prefix + "_filtered_locations_with_softclipping.bed"),(file_prefix + "_reduced_reference_genome.fa")]:
	if os.path.isfile(res_file):
		os.remove(res_file)

# filter BAM file for reads with a mapping quality lower than 20
bam_filt = file_prefix + "_filtered.bam"
subprocess.call("samtools view -b -h -q 20 " + usr_args.bam_file + " > " + bam_filt, shell=True)
subprocess.call("samtools index " + bam_filt, shell=True)

# create bed file from bam file, merge overlapping/touching/12bp apart read locations and filter on region length (set at 200bp)
subprocess.call("bedtools bamtobed -split -i " + bam_filt + " | bedtools merge -d 12 | awk \'{if (($3-$2)>80){print $0}}\' | sort -k1,1 -k2,2n > " + file_prefix + "_filtered.bed", shell=True)

# intersect filtered hits with known alpharetroviral locations in the genome (identified with the run_blast_ref_seq.sh script)
subprocess.call("bedtools intersect -v -a " + file_prefix + "_filtered.bed -b " + usr_args.ref_pos + " > " + file_prefix + "_filtered_locations.bed", shell=True)
#subprocess.call("rm " + file_prefix + "_filtered.bed", shell=True)

# extract all intervals of interest from BAM file and return used intervals to file
if os.stat(file_prefix + "_filtered_locations.bed").st_size > 0:
	sc_ints = open(file_prefix + "_filtered_locations_with_softclipping.bed", "a")
	for line in open(file_prefix + "_filtered_locations.bed").read().rstrip("\n").split("\n"):
		# process interval into searchable string
		interval = str(line.split("\t")[0]) + ":" + str(line.split("\t")[1]) + "-" + str(line.split("\t")[2])
		# for the reads within each interval, extract the CIGAR scores, remove numbers and duplicates, then store in a list
		cigar_scores = (subprocess.check_output("samtools view " + bam_filt + " " + interval + " | cut -f6 | tr -d \'0123456789\' | sort | uniq", shell=True)).rstrip("\n").split("\n")
		read_count = int((subprocess.check_output("samtools view " + bam_filt + " " + interval + " | wc -l", shell=True)).rstrip("\n"))
		# if there are cases of mapped followed by soft-clipped AND soft-clipped followed by mapped, extract the reads for that interval to file
		if read_count > 3:
			if all(x in cigar_scores for x in ["MS", "SM"]):
				subprocess.call("samtools view " + bam_filt + " " + interval + " >> " + file_prefix + "_filtered_locations.sam", shell=True)
				sc_ints.write(line + "\t1\t" + line_prefix + "\t+\n")
			elif any(x in cigar_scores for x in ["MS", "SM"]):
				subprocess.call("samtools view " + bam_filt + " " + interval + " >> " + file_prefix + "_filtered_locations.sam", shell=True)
				sc_ints.write(line + "\t2\t" + line_prefix + "\t+\n")
		else:
			subprocess.call("rm " + file_prefix + "_filtered_locations.bed " + file_prefix + "_filtered.bam*", shell=True)
			sys.exit(file_prefix.replace("_viral_mapped","") + " - No identified locations not in reference list.")
	sc_ints.close()

else:
	subprocess.call("rm " + file_prefix + "_filtered_locations.bed " + file_prefix + "_filtered.bam*", shell=True)
	sys.exit(file_prefix.replace("_viral_mapped","") + " - No identified locations not in reference list.")
	
