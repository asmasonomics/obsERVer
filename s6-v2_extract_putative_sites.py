## python s6-v2_extract_putative_sites.py [-h] viral_mapped.sorted.bam

import argparse
import os
import subprocess
from collections import Counter

parser = argparse.ArgumentParser(
	description="s6-v2_extract_putative_sites.py takes the viral_mapped bam file, converts to bed and then extracts putative insert regions",
	epilog="Author: Andrew Mason; Release: 26/10/2017; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("bam_file", help="specify viral mapped BAM file for analysis")
usr_args = parser.parse_args()

file_prefix = ((usr_args.bam_file).split("/")[-1]).split(".")[0]
line_prefix = (((usr_args.bam_file).split("_viral")[0]).split("_")[-1]).replace("line","").replace("indiv","")
#ref_pos = "/exports/cmvm/eddie/eb/groups/CTLGH_GCRF/ALVEs_ASM/galGal5_alpha_pos.bed"
ref_pos = "/exports/cmvm/eddie/eb/groups/CTLGH_GCRF/ALVEs_ASM/05_gal6_work/g6_alpha_pos.bed"
pseudo = "/exports/cmvm/eddie/eb/groups/CTLGH_GCRF/ALVEs_ASM/pseudochromosome_seq.fa"
alve_loc = "/exports/cmvm/eddie/eb/groups/CTLGH_GCRF/ALVEs_ASM/ALVE_locations.bed"


# check there are no previous run output files
for res_file in [(file_prefix + "_filtered_locations.sam"),(file_prefix + "_filtered_locations_with_softclipping.bed"),(file_prefix + "_reduced_reference_genome.fa")]:
	if os.path.isfile(res_file):
		os.remove(res_file)
		
# filter BAM file for reads with a mapping quality lower than 20
bam_filt = file_prefix + "_filtered.bam"
subprocess.call("samtools view -b -h -q 20 " + usr_args.bam_file + " > " + bam_filt, shell=True)
subprocess.call("samtools index " + bam_filt, shell=True)

# create bed file from bam file, merge overlapping/touching/12bp apart read locations and filter on region length (set at 150bp) **edit, to reduce false positives
subprocess.call("bedtools bamtobed -split -i " + bam_filt + " | bedtools merge -d 12 | awk \'{if (($3-$2)>150){print $0}}\' | sort -k1,1 -k2,2n > " + file_prefix + "_filtered.bed", shell=True)

# intersect filtered hits with known alpharetroviral locations in the genome (identified with the run_blast_ref_seq.sh script)
subprocess.call("bedtools intersect -v -a " + file_prefix + "_filtered.bed -b " + ref_pos + " > " + file_prefix + "_filtered_locations.bed", shell=True)
subprocess.call("rm " + file_prefix + "_filtered.bed", shell=True)

# extract all intervals of interest from BAM file and return used intervals to file
sc_ints = open(file_prefix + "_filtered_locations_with_softclipping.bed", "a")
for line in open(file_prefix + "_filtered_locations.bed").read().rstrip("\n").split("\n"):
	# process interval into searchable string
	interval = str(line.split("\t")[0]) + ":" + str(line.split("\t")[1]) + "-" + str(line.split("\t")[2])
	# for the reads within each interval, extract the CIGAR scores, remove numbers and duplicates, then store in a list
	cigar_scores = (subprocess.check_output("samtools view " + bam_filt + " " + interval + " | cut -f6 | tr -d \'0123456789\' | sort | uniq", shell=True)).rstrip("\n").split("\n")
	read_count = int((subprocess.check_output("samtools view " + bam_filt + " " + interval + " | wc -l", shell=True)).rstrip("\n"))
	
	# empty variables used for file output and internal modifications
	hex_location = ""
	hom_res = ""
	softclip_sup = ""
	hex_st = []
	hex_end = []
	ALVE_overlap = ""
	
	if (any(x in cigar_scores for x in ["MS", "SM"]) and (read_count > 1)):
		# generate interval and put in reduced sam file
		if (int((subprocess.check_output("samtools view " + bam_filt + " " + interval + " | wc -l ", shell=True)).rstrip("\n")) < 2):
			break
		subprocess.call("samtools view " + bam_filt + " " + interval + " >> " + file_prefix + "_filtered_locations.sam", shell=True)
		
		# generate temporary files based on cigar strings for reads supporting the 5' and 3' ends of an insertion
		subprocess.call("samtools view " + bam_filt + " " + interval + " | awk \'$6~/^[0-9]+M[0-9]+S$/{split($6,a,\"M\"); print $3 \"\t\" $4 \"\t\" a[1]\"M\t\" a[2] \"\t\" $10}\' > " + file_prefix + "_tmpM", shell=True)
		subprocess.call("samtools view " + bam_filt + " " + interval + " | awk \'$6~/^[0-9]+S[0-9]+M$/{split($6,a,\"S\"); print $3 \"\t\" $4 \"\t\" a[1]\"S\t\" a[2] \"\t\" $10}\' > " + file_prefix + "_tmpS", shell=True)
		
		
		# write interval to file in BED format and overlap with list of known ALVEs, and extract information for outfile
		check_int = open(file_prefix + "_check_int.bed","w")
		check_int.write(interval.split(":")[0] + "\t" + str((interval.split(":")[1]).split("-")[0]) + "\t" + str((interval.split(":")[1]).split("-")[1]))
		check_int.close()
		int_overlap = ((subprocess.check_output("bedtools intersect -wo -a " + file_prefix + "_check_int.bed -b " + alve_loc, shell=True)).rstrip("\n")).split("\n")
		
		# if there are any overlapping sequences
		if int_overlap[0]:
			# if there is 1 
			if len(int_overlap) == 1:
				ALVE_overlap = str((int_overlap[0]).split("\t")[6]) + "-" + (str((int_overlap[0]).split("\t")[4])) + "-" + (str((int_overlap[0]).split("\t")[5]))
			# if there are many create a ; separated list of possible overlaps
			if len(int_overlap) > 1:
				for p in int_overlap:
					ALVE_overlap += str((int_overlap[0]).split("\t")[6]) + "-" + (str((int_overlap[0]).split("\t")[4])) + "-" + (str((int_overlap[0]).split("\t")[5])) + ";"
				ALVE_overlap = ALVE_overlap[:-1]
		else:
			ALVE_overlap = "."
						
		# check sequences which give 5prime support for the ALVE insertion
		if os.stat(file_prefix + "_tmpM").st_size != 0:
			softclip_sup = "5pri"
			SC_len = []
			# for each read
			for i in open(file_prefix + "_tmpM").read().rstrip("\n").split("\n"):
				# add the softclipped length to list
				SC_len += [(int((i.split("\t")[3]).replace("S","")))]
				# if the soft clipped length is greater than 12bp write the sequence to file and perform a short blastn against the pseudo to check for ALVE homology
				if (int((i.split("\t")[3]).replace("S",""))) > 12:
					check_seq = (i.split("\t")[-1])[(int((i.split("\t")[2]).replace("M",""))-1):]
					seq = open(file_prefix + "_check_seq.fa","w")
					seq.write(">check\n" + check_seq)
					seq.close()
					if (subprocess.check_output("blastn -db " + pseudo + " -query " + file_prefix + "_check_seq.fa -task \"blastn-short\" -outfmt \"6\" -evalue 0.0001 | wc -l",shell=True).rstrip("\n")) > 0:
						hex_end += [(int(i.split("\t")[1])+(int((i.split("\t")[2]).replace("M",""))-1))]
			
			# if all SC lengths are shorter than 13bp designate the short length unsuitable for BLAST check
			SC_len.sort()
			if (len(SC_len) == (sum(k < 13 for k in SC_len))):
				hom_res = "5priSHORT"
			if (len(hex_end) > 0):
				hom_res = "5priYES"
			
					
				
		# check those sequences which give 3prime end support, as above		
		if os.stat(file_prefix + "_tmpS").st_size != 0:		
			if not softclip_sup:
				softclip_sup = "3pri"
			else:
				softclip_sup += "-3pri"
			SC_len = []
			for j in open(file_prefix + "_tmpS").read().rstrip("\n").split("\n"):
				SC_len += [(int((j.split("\t")[2]).replace("S","")))]
				if (int((j.split("\t")[2]).replace("S",""))) > 12:
					check_seq = (j.split("\t")[-1])[0:(int((j.split("\t")[2]).replace("S","")))]
					seq = open(file_prefix + "_check_seq.fa","w")
					seq.write(">check\n" + check_seq)
					seq.close()
					if (subprocess.check_output("blastn -db " + pseudo + " -query " + file_prefix + "_check_seq.fa -task \"blastn-short\" -outfmt \"6\" -evalue 0.0001 | wc -l",shell=True).rstrip("\n")) > 0:
						hex_st += [int(j.split("\t")[1])]
				
			SC_len.sort()
			if (len(SC_len) == (sum(k < 13 for k in SC_len))):
				if not hom_res:
					hom_res = "3priSHORT"
				else:
					hom_res += "-3priSHORT"
			if (len(hex_st) > 0):
				hom_res += "-3priYES"
				
			
		# remove duplicates from list	
		st2 = list(set(hex_st))
		st2.sort()
		st3 = list(set(hex_end))
		st3.sort()
		
		# if there is one option for 5pri and 3pri, check they satisfy hexamer spacing 
		if ((len(st2) == 1) and (len(st3) == 1)):
			if ((st2[0] < st3[0]) and ((st2[0]+10) > st3[0])):
				hex_location = str(st2[0]) + "-" + str(st3[0])
			elif (st2[0] == st3[0]):
				hex_location = "?-" + str(st2[0]) + "-?"
			else:
				if (len(hex_st) > len(hex_end)):
					hex_location = str(st2[0]) + "-?"
				elif (len(hex_st) < len(hex_end)):
					hex_location = "?-" + str(st3[0])
				else:
					hex_location = "."
		
		# if both positions lists have one or multiple options iterate through to find pairs that satisfy the hexamer sizings.
		elif ((len(st2) >= 1) and (len(st3) >= 1)):
			if ((len(st2) == 1) and (len(st3) > 1)):
				for k in st3:
					if ((k >= st2[0]) and (k < (st2[0] + 10))):
						hex_location = str(st2[0]) + "-" + str(k)
				if hex_location == "":
					most_common = [k for k in Counter(hex_end).most_common()]
					if (int(most_common[0][1]) > int(most_common[1][1])):
						if (int(most_common[0][1]) > len(hex_st)):
							hex_location = "?-" + str(most_common[0][0])
						elif (int(most_common[0][1]) < len(hex_st)):
							hex_location = str(st2[0]) + "-?"
						else:
							hex_location = "."
					elif (int(most_common[0][1]) == int(most_common[1][1])):
						if (int(most_common[0][1]) > len(hex_st)):
							hex_location = "?-" + str(most_common[0][0])
						elif (int(most_common[0][1]) < len(hex_st)):
							hex_location = str(st2[0]) + "-?"
						else:
							hex_location = "."
					else:
						hex_location = "." 
			
			elif ((len(st3) == 1) and (len(st2) > 1)):
				for k in st2:
					if ((k <= st3[0]) and (k > (st3[0] - 10))):
						hex_location = str(k) + "-" + str(st3[0])
				if hex_location == "":
					most_common = [k for k in Counter(hex_st).most_common()]
					if (int(most_common[0][1]) > int(most_common[1][1])):
						if (int(most_common[0][1]) > len(hex_end)):
							hex_location = str(most_common[0][0]) + "-?"
						elif (int(most_common[0][1]) < len(hex_end)):
							hex_location = "?-" + str(st3[0])
						else:
							hex_location = "."
					elif (int(most_common[0][1]) == int(most_common[1][1])):
						if (int(most_common[0][1]) > len(hex_end)):
							hex_location = str(most_common[0][0]) + "-?"
						elif (int(most_common[0][1]) < len(hex_end)):
							hex_location = "?-" + str(st2[0])
						else:
							hex_location = "."
					else:
						hex_location = "."
					
			elif ((len(st2) > 1) and (len(st3) > 1)):
				for p in st2:
					for q in st3:
						if ((p < q) and ((p + 10) >= q)):
							hex_location = str(p) + "-" + str(q)
				
				# if there is no hex satisfying pair, attempt to give most supported end using multiple occurences
				if hex_location == "":
					st2_most_common = [k for k in Counter(hex_st).most_common()]
					st3_most_common = [k for k in Counter(hex_end).most_common()]
					
					if ((int(st2_most_common[0][1]) > int(st2_most_common[1][1])) and (int(st3_most_common[0][1]) > int(st3_most_common[1][1]))):
						if (int(st2_most_common[0][1]) > int(st3_most_common[0][1])):
							hex_location = str(st2_most_common[0][0]) + "-?"
						elif (int(st2_most_common[0][1]) < int(st3_most_common[0][1])):
							hex_location = "?-" + str(st3_most_common[0][0])
						else:
							hex_location = "."
					
					elif ((int(st2_most_common[0][1]) == int(st2_most_common[1][1])) and (int(st3_most_common[0][1]) == int(st3_most_common[1][1]))):
						hex_location = "."
					
					elif ((int(st2_most_common[0][1]) > int(st2_most_common[1][1])) and (int(st3_most_common[0][1]) == int(st3_most_common[1][1]))):
						if (int(st2_most_common[0][1]) > int(st3_most_common[0][1])):
							hex_location = str(st2_most_common[0][0]) + "-?"
						else:
							hex_location = "."
					
					elif ((int(st2_most_common[0][1]) == int(st2_most_common[1][1])) and (int(st3_most_common[0][1]) > int(st3_most_common[1][1]))):
						if (int(st3_most_common[0][1]) > int(st2_most_common[0][1])):
							hex_location = "?-" + str(st3_most_common[0][0])
						else:
							hex_location = "."
					
					else:
						hex_location = "."
		
		# if there is one possibility for one side and the other is empty
		elif ((len(st2) == 1) and (len(st3) == 0)):
			hex_location = str(st2[0]) + "-?"					
						
		# as above, but alternate
		elif ((len(st2) == 0) and (len(st3) == 1)):
			hex_location = "?-" + str(st3[0])
		
		# if there are no designations for 5pri, but muliple possibilities for 3pri, check for most common 3pri and use that. If there are equal most common designate as null.
		elif ((len(st2) == 0) and (len(st3) > 1)):
			most_common = [k for k in Counter(hex_end).most_common()]
			if (int(most_common[0][1]) > int(most_common[1][1])):
				hex_location = "?-" + str(most_common[0][0])
			else:
				hex_location = "."
		
		# as above but alternate scenario
		elif ((len(st2) > 1) and (len(st3) == 0)):
			most_common = [k for k in Counter(hex_st).most_common()]
			if (int(most_common[0][1]) > int(most_common[1][1])):
				hex_location = str(most_common[0][0]) + "-?"
			else:
				hex_location = "."
		
		
		# fill in gap variables if no designations have been made
		elif ((len(st2) == 0) and (len(st3) == 0)):
			hex_location = "." 
		
		if ((os.stat(file_prefix + "_tmpS").st_size == 0) and (os.stat(file_prefix + "_tmpM").st_size == 0)):
			softclip_sup == "."
		if not hom_res:
			hom_res = "."
		if not hex_location:
			hex_location = "."
		
		# write final output file
		# File is in BED12 format with the fifth, eleventh and twelfth columns blank (".")
		# C1 = contig, C2 = start, C3 = end, C4 = data prefix, C6 = strand ("+" in all cases), C7 = known ALVE overlap, C8 = softclip read support, C9 = ALVE BLAST homology result, C10 = putative hex location (st-end)
		# C7 gives the name of any overlapped feature, its hexamer and hexamer location - matches are done by full region overlap, so could include multiple overlaps (given as a ; separated list)
		# C8 states whether there are softclipped reads in the given interval, and whether there are 5pri, 3pri or both (5pri-3pri)
		# C9 gives an indication from the ALVE BLAST against the pseudo. This will be empty (".") if no reads gave support, but were long enough for a BLAST. If there were no reads longer than 12bp (for good BLAST matches) for either 5pri or 3pri it will show 5priSHORT. Good matches will have 5priYES or similar for the 3pri
		# C10 gives a predicted hexamer location from the clipped reads - this is based on the best supported sites and whether there are both 5pri  and 3pri designations. If 5pri is confident but 3pri missing the 3pri will be represented as "-?" etc.
		# All these designations give extra info, but does not replace manual checking. Also. Hard clipped reads are not included, but should be checked, especially if there is too short sequence for ALVE BLAST or if designing assays.
		sc_ints.write(line + "\t" + line_prefix + "\t.\t+\t" + ALVE_overlap + "\t" + softclip_sup + "\t" + hom_res.replace("--","-") + "\t" + hex_location + "\t.\t.\n")
	

sc_ints.close()

#remove extra files and clear up
rem_files = [bam_filt, (bam_filt + ".bai"), (file_prefix + "_filtered_locations.bed"), (file_prefix + "_tmpM"), (file_prefix + "_tmpS"), (file_prefix + "_check_seq.fa"), (file_prefix + "_check_int.bed")]
for z in rem_files:
	if os.path.exists(z):
		os.remove(z)
