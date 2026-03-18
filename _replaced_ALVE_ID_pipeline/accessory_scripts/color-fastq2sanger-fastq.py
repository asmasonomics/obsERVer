## python color-fastq2sanger-fastq.py [-h] file.fastq.gz

import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(
        description="convert SOLiD color fastq to 'standard' sanger fastq (phred33 is existing qual score)",
        epilog="Author: Andrew Mason; Release: 06/11/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("fastq", help="fastq file for conversion (.fastq.gz format)")
usr_args = parser.parse_args()

fastq_out = (usr_args.fastq).replace(".gz","").replace(".fastq","_converted2sanger.fastq")

color_dict = {'T0':'T', 'T1':'G', 'T2':'C', 'T3':'A', 'T.':'N', 'C0':'C', 'C1':'A', 'C2':'T', 'C3':'G', 'C.':'N', 'G0':'G', 'G1':'T', 'G2':'A', 'G3':'C', 'G.':'N', 'A0':'A', 'A1':'C', 'A2':'G', 'A3':'T', 'A.':'N', 'N0':'N', 'N1':'N', 'N2':'N', 'N3':'N', 'N.':'N'}

fastq = []
tmp = []
i = 0
seq = ""

for line in ((subprocess.check_output("zcat " + usr_args.fastq, shell=True)).rstrip("\n")).split("\n"):
	i += 1
	
	if i==2:
		seq += color_dict[line[0:2]]
		for val in line[2:]:
			seq += color_dict[(seq[-1] + val)]	
		tmp += [seq]
		seq = ""
	elif i==4:
		tmp += [line[1:]]
		fastq += [tmp]
		tmp = []
		i = 0
	else:
		tmp += [line]

out = open(fastq_out, "w")
for i in fastq:
	for j in i:
		out.write(j + "\n")
out.close()
subprocess.call("gzip " + fastq_out, shell=True)

