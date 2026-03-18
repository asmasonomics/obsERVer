## python seq_extract.py [-h] [--prefix PREFIX] positions_file ref_genome

import argparse
parser = argparse.ArgumentParser(
	description="seq_extract extracts sequences from a reference fasta file using given positional information.",
	epilog="Author: Andrew Mason; Release: 04/11/14; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("pos_file", help="Positions file to be analysed - column order (tab spaced): chr st end str source_program")
parser.add_argument("ref_file", help="Reference fasta file, often a reference genome")
parser.add_argument("--prefix", default="USER", help="Specify your own output file prefix (default = \"USER\")")
parser.add_argument("-u", "--use_full", action="store_true", help="if select seq-extract will use full contig name for extract, rather than default not")
usr_args = parser.parse_args()

import sys
import static_functions

# Create dict of reference genome for sequence extraction
ref_dict = static_functions.seq_dict_creator(usr_args.ref_file)

# Open positions file and format list
pos_list = static_functions.list_initial_formatter(usr_args.pos_file)

#create a file for extracted sequences to be written to and extract sequences
out_name = "./" + usr_args.prefix + "_seq.fasta"
cont_full = "N"
if usr_args.use_full:
    cont_full = "Y"
static_functions.seq_extract(out_name, pos_list, ref_dict, cont_full)

print("Sequence extract complete.\n")
