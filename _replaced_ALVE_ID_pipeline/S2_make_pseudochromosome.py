## python make_pseudochromosome.py [-h] [--gaps GAP(INT)] [--ref_gen REF_GENOME.FA] ref_seq.fa

import argparse
import subprocess

parser = argparse.ArgumentParser(
	description="make_pseudochromosome takes a reference fasta file and creates a single pseudochromosome it then add to a given genome.",
	epilog="Author: Andrew Mason; Release: 17/11/15; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ref_seq", help="fasta of reference sequences used to create the pseudochromosome")
parser.add_argument("--ref_gen", help="reference genome to have pseudochromosome added to")
parser.add_argument("--gaps", type=int, default=500, help="<INT> Specify gap in bases between reference sequence on pseudochromosome (default = 500)")
usr_args = parser.parse_args()

# ensure input fasta file has sequences on a single line
subprocess.call(("awk \'/^>/ {print (NR==1 ? \"\" :RS)$0; next}{printf \"%s\", $0}END{printf RS}\' " + usr_args.ref_seq + " > ref_seq.tmp"), shell=True)

# define lists and separate reference sequence file into headers and sequences
headers = []
seq = []
lengths = []

for line in open("ref_seq.tmp").read().rstrip("\n").split("\n"):
	if line.startswith(">"):
		headers += [line.replace(">","")]
	else:
		seq += [line]
		lengths += [len(line)]

# prepare pseudochromosome reference sequence positions file and define initial variables
out_pos = open("pseudochromosome_positions.txt", "w")
# pseudochromosome starts with specified gap (default = 500)
pseudo = ("N" * usr_args.gaps)
pos = usr_args.gaps + 1
i = 0
for header in headers:
	out_pos.write(header + "\t" + str(pos) + "\t" + str((pos - 1) + lengths[i]) + "\t" + str(lengths[i]) + "\n")
	pos += (lengths[i] + usr_args.gaps)
	# each reference sequence is separated by specified gap (default = 500)
	pseudo += seq[i] + ("N" * usr_args.gaps)
	i += 1

# pseudochromosome also finishes with specified gap (default = 500)
pseudo += ("N" * usr_args.gaps)
out_pos.close()
# write pseudochromosome sequence to file
out_seq = open("pseudochromosome_seq.fa", "w")
out_seq.write(">pseudochromosome 1:" + str(len(pseudo)) + "\n" + pseudo + "\n")
out_seq.close()

if usr_args.ref_gen:
	# concatenate pseudochromosome with reference genome and index for mapping (long last step)
	subprocess.call(("awk \'/^>/ {print (NR==1 ? \"\" :RS)$0; next}{printf \"%s\", $0}END{printf RS}\' " + usr_args.ref_gen + " > ref_gen.tmp"), shell=True)
	outname = (((usr_args.ref_gen).split(".")[0]).split("/")[-1]) + "_plus_pseudochromosome.fa"
	subprocess.call("cat pseudochromosome_seq.fa ref_gen.tmp > " + outname + " ; rm ref_seq.tmp ref_gen.tmp", shell=True)
	subprocess.call("bwa index " + outname, shell=True)
else:
	subprocess.call("bwa index pseudochromosome_seq.fa; rm ref_seq.tmp", shell=True)


