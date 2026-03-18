## static_functions.py - a module of regularly used functions

import re
import os
import sys
import subprocess




## FUNCTION 1 - takes a multi-fasta file and returns a list of all the headers ##
def header_extractor(fasta_file):
        headers = []
        for line in (fasta_file.split("\n")):
                if line.startswith(">"):
                        headers += [(re.sub(r"[ ]+", "_", (line.replace(">", ""))))]
        return headers


## FUNCTION 2 - takes a multi-fasta file, extracts the sequences, rather than headers, and returns an ordered list of sequences ##
def seq_only_extractor(fasta_file):
        seq_list = re.split(r"_", (re.sub(r"[\r|\n]", "", (re.sub(r"(>(.+)\n)", "_", fasta_file)))))
        seq_list.pop(0)
        return seq_list


## FUNCTION 3 - creates a callable dictionary for extracting sequences ##
def seq_dict_creator(ref_fasta_file):
        ref_file = open(ref_fasta_file).read().rstrip("\n")
        ref_headers = header_extractor(ref_file)
        ref_seq_list = seq_only_extractor(ref_file)
        ref_dict = {}
        i=0
        for name in ref_headers:
                        ref_dict[name] = ref_seq_list[i]
                        i+=1
                        
        return ref_dict


## FUNCTION 4 - extracts sequences from reference genome given a set of positions in the required format ##
def seq_extract(output_file_name, positions_list, genome_dict, yn):
        seq_file = open(output_file_name, "w")
        for pos in positions_list:
                chr_name = str(pos[0]) + "_"
                if yn=="Y":
                        chr_name = (str(pos[0])).replace("_1", ".1")
                sequence = next(val for key, val in genome_dict.iteritems() if key.startswith(chr_name))
                if pos[2] > pos[1]:
                        start = pos[1] - 1
                        end = pos[2]
                else:
                        start = pos[2] - 1
                        end = pos[1]
		
		if yn=="Y":
                        chr_name += "_"
		
                if str(pos[3]) in ('+', 'm'):
                        hit_seq = sequence[start:end]
                        source = pos[4]
                        seq_file.write(">" + chr_name + str(start + 1) + ":" + str(end) + "_" + (pos[3]) + "_" + source + "\n" + hit_seq + "\n")
                else:
                        hit_seq = sequence[start:end]
                        source = pos[4]
                        lower_rev = hit_seq.replace("A" , "t").replace("T" , "a").replace("G" , "c").replace("C" , "g")
                        rev_comp = lower_rev.upper()
                        reversed_comp = rev_comp[::-1]
                        seq_file.write(">" + chr_name + str(start + 1) + ":" + str(end) + "_" + (pos[3]) + "_" + source + "\n" + reversed_comp + "\n")	
        seq_file.close() 
        return

## FUNCTION 5 - split multifasta into individual query files ##
def fasta_splitter(fasta_file, colon_ok):
        headers = header_extractor(fasta_file)
        seq_list = seq_only_extractor(fasta_file)
        i=0
        query_seq_files = []
        for seq in seq_list:
                seq_name = (headers[i] + ".fas")
                if colon_ok == "N":
                        seq_name = (headers[i] + ".fas").replace(":", "-")
                query_seq_files += [seq_name]
                query = open((seq_name), "w")
                query.write(">" + (headers[i]) + "\n" + seq + "\n")
                query.close()
                i+=1
        return query_seq_files


## FUNCTION 5a - split multifasta into individual query files (FASTER) ##
## use if you have already generated headers and seq_list in your main code ##
## does not return list of headers as above in function 5 ##
def fast_fasta_splitter(headers, seq_list, colon_ok):
        i=0
        for seq in seq_list:
                seq_name = (headers[i] + ".fas")
                if colon_ok == "N":
                        seq_name = (headers[i] + ".fas").replace(":", "-")
                query = open((seq_name), "w")
                query.write(">" + (headers[i]) + "\n" + seq + "\n")
                query.close()
                i+=1
        return 
                

## FUNCTION 6 - takes a sorted positions file and processes into an iterable list ##
def list_formatter(sorted_pos, sort_by):
	positions_split = sorted_pos.rstrip("\n").split("\n")
	pos_list = []
	by_chro = []
	#this loop splits up the list to create sublists and alters pos1 + pos2 to int()
	i=0
	for element in positions_split:
		x = (element.split("\t"))
		positions = []
                # this if block gets the start and end in the same orientation irrespective of strand +/-
		if int(x[1]) > int(x[2]):
			st, end = 1, 2
			x[end], x[st] = x[st], x[end]
		count = 0
		for value in x:
			count += 1
			if count == 2 or count == 3:
				y = int(value)
				positions += [y]
			#not changing to int for pos[0] or pos[3] which may use non number characters
			else:
				positions += [value]
		
		if i==0: 
			by_chro += [positions]
		elif positions[sort_by] == (by_chro[-1])[sort_by]:
			by_chro += [positions]
		else:
			pos_list += [by_chro]
			by_chro = []
			by_chro += [positions]
			
		if i == (len(positions_split) - 1):
			pos_list += [by_chro]
			
		i+=1
	return pos_list



## FUNCTION 7 - reduced function for list_formatter ##
def list_initial_formatter(pos_file):
		data = ((open(pos_file).read().rstrip("\n")).replace("\t", ",")).split("\n")
                pos_list = []
		for element in data:
                        x = (element.split(","))
			positions = []
			# this if block gets the start and end in the same orientation irrespective of strand +/-
			if int(x[1]) > int(x[2]):
				st, end = 1, 2
				x[end], x[st] = x[st], x[end]
			count = 0
			for value in x:
				count += 1
				if count == 2 or count == 3:
					y = int(value)
					positions += [y]
				#not changing to int for pos[0] or pos[3] which may use non number characters
				else:
					positions += [(value.replace("_1", ".1"))]
			pos_list += [positions]
		return pos_list


## FUNCTION 8 - takes a list of positions and merges them according to set requirements ##
def positions_merger(x, log_file):	
	x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))
	i = 0
	j = 0
	k = 0
	sort = None
	
	contig = x[0][0]
	contig_pos_number = str(len(x))
	#print("Processing positions from contig \"" + contig + "\". " + str(len(x)) + " position(s) to process.\n")
	#log_file.write("Processing positions from contig \"" + contig + "\". " + str(len(x)) + " position(s) to process.\n")
	
	while i < len(x) - 1:
		if sort == True:
			x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))
	
		# same contig/chromosome, i+1 end pos is within the i element --> remove i+1
		if (x[i][0] == x[i + 1][0]) and (x[i + 1][2] <= x[i][2]):
			x.pop(i+1)
			k+=1
			sort=True 
		
		# same contig/chromosome and i end pos is greater than i+1 start pos OR elements touch OR there is less than 11bp observed gap between elements 
		elif (x[i][0] == x[i + 1][0]) and ((x[i][2] >= x[i + 1][1]) or (x[i + 1][1] == ((x[i][2]) + 1)) or ((int(x[i+1][1]) - int(x[i][2])) <= 10)):
			#merged_pos = str(x[i][0]) + ", " + str(x[i][1]) + ", " + str(x[i+1][2]) + ", m, " +  str(x[i][4]) + "-" + str(x[i+1][4])
                        merged_pos = str(x[i][0]) + ", " + str(x[i][1]) + ", " + str(x[i+1][2]) + ", " + str(x[i][3]) + str(x[i+1][3]) + ", " + str(x[i][4]) + "," + str(x[i+1][4]) 
			mp = merged_pos.split(", ")
			new_pos = []
			count = 0
			for element in mp:
				count+=1
				if count == 2 or count == 3:
					y = int(element)
					new_pos += [y]
				else:
					new_pos += [element]
			x += [new_pos]
			x.pop(i+1)
			x.pop(i)
			k+=1
			sort = True
		else:
			i+=1
			sort = False
	
		if (i%50==0) and (i>0) and (i>j):
			j = i
			#print("The first " + str(i) + " list positions have been processed.\n\tOriginal list length: " + contig_pos_number + "\n\t" + str(k) + " list items merged or removed\n\tCurrent list length: " + str(len(x)) + "\n")
			#log_file.write("The first " + str(i) + " list positions have been processed.\n\tOriginal list length: " + contig_pos_number + "\n\t" + str(k) + " list items merged or removed\n\tCurrent list length: " + str(len(x)) + "\n")

        for pos in x:
                # sort strands
                m_count = (str(pos[3])).count("m")
                plus_count = (str(pos[3])).count("+")
                minus_count = (str(pos[3])).count("-")
                strand_length = m_count + plus_count + minus_count
                new_strand = ""
                if strand_length > 1:
                        if plus_count >= minus_count:
                                new_strand = "+"
                        else:
                                new_strand = "-"
                        pos[3] = new_strand

                # sort source
                source_line = list(set((pos[4]).replace("-", ",").split(",")))
                source_line.sort()
                pos[4] = ("-".join(source_line)).rstrip(" ")
                
	return(x)


## FUNCTION 9 - processes sequences with RM annotation for CR1 elements and removes them where necessary
def cr1_fragment_remover(pre_pro_CR1_file, updated_pos_file, update_num, log_file, prefix):
        print("Removing CR1 annotated fragments.");log_file.write(" Removing CR1 annotated fragments. ")
        raw_CR1_data = (open(pre_pro_CR1_file).read().rstrip("\n")).split("\n\n")
        CR1_data = []
        # split data so each element of each individual fragment is available
        for a in raw_CR1_data:
                x = []
                for b in (a.split("\n")):
                        y = []
                        for c in (b.split("\t")):
                                y += [c]
                        x += [y]
                CR1_data += [x]

        new_hits = []
        for hit in CR1_data:
                if len(hit) > 1:
                        x = len(hit) - 1
                        #define key parts of fragment for rewriting
                        if ((hit[0][1].split(":"))[0]).count("_") > 1:
                                chromo = str(((hit[0][1]).split("_"))[0]) + "_" + str(((hit[0][1]).split("_"))[1])
                                start_pos = (((hit[0][1]).split("_"))[2].split(":"))[0]
                                end_pos = (((hit[0][1]).split("_"))[2].split(":"))[1]
                                strand = ((hit[0][1]).split("_"))[3]
                        else:
                                chromo = ((hit[0][1]).split("_"))[0]
                                start_pos = (((hit[0][1]).split("_"))[1].split(":"))[0]
                                end_pos = (((hit[0][1]).split("_"))[1].split(":"))[1]
                                strand = ((hit[0][1]).split("_"))[2]
                        #add original fragment in correct form so it can be removed easily using uniq -u
                        new_hits += [(str(chromo) + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + strand + "\t" + prefix)]
	
                        #processing when CR1 is the first annotated part of the fragment
                        if hit[0][5] == "LINE/CR1":
                                end = hit[0][3]
                                sw = hit[0][0]
                                # if the end of the CR1 annotated fragment overlaps the start of the LTR annotated section see which SW score is higher and give that annotation preference
                                if end > hit[1][2]:
                                        if int(sw) > int(hit[1][0]):
                                                new_start_pos = int(start_pos) + int(end)
                                                new_hits += [(str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\t" + prefix)]
			
                                        else:
                                                new_start_pos = int(start_pos) + int(hit[1][2])
                                                new_hits += [(str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\t" + prefix)]
                                else:
                                        new_start_pos = int(start_pos) + int(end)
                                        new_hits += [(str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\t" + prefix)]
	
                        #processing when CR1 is the last annotated part of the fragment
                        if hit[x][5] == "LINE/CR1":
                                st = hit[x][2]
                                sw = hit[x][0]
                                # if the start of the CR1 annotated fragment overlaps the end of the LTR annotated section see which SW score is higher and give that annotation preference
                                if st < hit[x-1][3]:
                                        if int(sw) > int(hit[x-1][0]):
                                                new_end_pos = int(start_pos) + (int(st) - 1)
                                                new_hits += [(str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\t" + prefix)]
			
                                        else:
                                                new_end_pos = int(start_pos) + int(hit[x-1][3])
                                                new_hits += [(str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\t" + prefix)]
                                else:
                                        new_end_pos = int(start_pos) + (int(st) - 1)
                                        new_hits += [(str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\t" + prefix)]

        print("CR1 fragments removed from putative sequence. Writing modified positions to file.");log_file.write("CR1 fragments removed from putative sequence. Writing modified positions to file.")
			
        #write newly edited hits to tmp file
        added_hits = open("CR1_removed.tmp", "w")
        for hit in new_hits:
                added_hits.write(hit + "\n")
        added_hits.close()

        # concatenate original positions list with edited positions, then sort the file and remove the duplicated positions (the original and pre-edit positions with CR1 annotation), then remove .tmp files
        subprocess.call("cat " + updated_pos_file + " CR1_removed.tmp > merge_plus_CR1_removed.tmp", shell=True)
        subprocess.call("sort -k1,1 -k2,2n merge_plus_CR1_removed.tmp | uniq -u > validated_positions_update_" + update_num + ".txt", shell=True)
        subprocess.call("rm *.tmp", shell=True)
        print("Modified positions successfully written to file.\n");log_file.write(" Modified positions successfully written to file.\n")

        return


## FUNCTION 10 - processes HMM files in folder into HMM flatfile databases
def run_hmm_press(hmm_location, database_root, log_file):
        if not os.path.isfile(hmm_location + database_root + ".h3i"):
                hmm_files = ""
                cwd = os.getcwd()
                os.chdir(hmm_location)
                for filename in os.listdir("./"):
                       if filename.endswith(".hmm"):
                                hmm_files += filename + " "
                subprocess.call("cat " + hmm_files + "> " + database_root, shell=True)
                subprocess.call("hmmpress " + database_root, shell=True)
                os.chdir(cwd)

                if log_file != "N":
                        print("Changing into defined HMM containing directory.");log_file.write("Changing into defined HMM containing directory.")
                        print("Concatenating HMM files to create \"" + database_root + "\" database flatfile.");log_file.write("Concatenating HMM files to create \"" + database_root + "\" database flatfile.")
                        print("Initiating hmmpress protocol to create binary indexed HMM database.");log_file.write("Initiating hmmpress protocol to create binary indexed HMM database.")
                        print("Successfully created \"" + database_root + "\" index files. Returning to original working directory.\n");log_file.write("Successfully created \"" + database_root + "\" index files. Returning to original working directory.\n")
        else:
                if log_file != "N":
                        print("\"" + database_root + "\" database already indexed. Proceeding with validation.\n");log_file.write("\"" + database_root + "\" database already indexed. Proceeding with validation.\n")

        return

















