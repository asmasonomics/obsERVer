### python seq_formatter.py [-h] [-rc REV_COMP] input_seq.fasta

import argparse
parser = argparse.ArgumentParser(
	description="seq_formatter is a multiple function script for splitting multi fasta files. It also produces reverse complement sequences (user option). If you have a single sequence fasta file it will detect this and automatically do a reverse complement.",
	epilog="Author: Andrew Mason; Release: 04/11/14; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("input_seq", help="Fasta file to be analysed")
parser.add_argument("-rc", "--rev_comp", help="Specify this option to also reverse complement all sequences in the specified input fasta, producing individual rc sequence files", action="store_true")
usr_args = parser.parse_args()

import sys
import re
import subprocess
import static_functions

fasta_format_checker = int(subprocess.check_output("grep \"^>\" " + usr_args.input_seq + " | wc -l | awk \'{print $1}\'", shell=True))
if (fasta_format_checker == 0):
        sys.exit("Input file not in fasta format. Ensure sequences have headers starting with \">\".")
        
# Open user input file, 
seq_file = open(usr_args.input_seq).read().rstrip("\n")
headers = static_functions.header_extractor(seq_file)
seq_list = static_functions.seq_only_extractor(seq_file)

# Check the default processing based on number of sequences in the user input file
splitter = True
if (len(headers) == 1):
        print("One sequence detected. Creating a reverse complement sequence.")
        usr_args.rev_comp = True
        splitter = False
        
# Run the fasta splitter function if required
print("Creating individual fasta files for each sequence.")
if (splitter == True):
        static_functions.fast_fasta_splitter(headers, seq_list, "Y")

# Run reverse complement if required
if (usr_args.rev_comp == True):
        print("Creating files of the reverse complement of the individual sequences.")
        rc_seq_list = []
        for element in seq_list:
                lower_rev = element.replace("A" , "t").replace("T" , "a").replace("G" , "c").replace("C" , "g")
                rev_comp = (lower_rev.upper())[::-1]
                rc_seq_list += [rev_comp]

        name = (re.sub(r"\.(.+)", "", sequence_file))
        file_name = name + "_rc.fasta"

        out_file = open(file_name, "w")
        i=0
        counter = (len(headers) - 1)
        while i < counter:
                out_file.write((headers[count]) + "\n" + rc[i] + "\n")
                i+=1
        out_file.close

print("Process complete.\n")
