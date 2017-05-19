from Bio.Blast import NCBIXML
import argparse
import re
from collections import defaultdict
import sys
import natsort
import os

parser = argparse.ArgumentParser(description='Take blastx output and detects linked mutations from each hit')
parser.add_argument('-f', '--fasta', help='FASTA file containing raw reads') #, required=True
parser.add_argument('-x', '--xml', help='BLASTx XML output file')
parser.add_argument('-e', '--evalue', help='E value cut off to filter out hits above this threshold', default=0.01)
parser.add_argument('-l', '--linked', action='store_true', help='flag to set whether linked output is enabled')
# parser.add_argument('-c', '--linkedcount', type=int, help='If linked output is enabled then only reads strands with this number of linked mutations or more will be outputted', default=1)
parser.add_argument('-u', '--unk', action='store_true', help='flag to set whether unknown mutations are outputted')
parser.add_argument('-p', '--pi', help='Protease mutation dictionary file, defaults to included file (PI_mutations.txt) if not specified', default="PI_mutations.txt")
parser.add_argument('-n', '--nrti', help='NRTI mutation dictionary file, defaults to included file (NRTI_mutations.txt) if not specified', default="NRTI_mutations.txt")
parser.add_argument('-nn', '--nnrti', help='NNRTI mutation dictionary file, defaults to included file (NNRTI_mutations.txt) if not specified', default="NNRTI_mutations.txt")

args = parser.parse_args()

#TODO: Include a paired end mode (or do in separate script after) that will combine mutations and counts for paired reads. Best way to do this would probably be to parse the BLASTx xml and store all the hits in a dictionary for each read then loop through the dictionary and find the matching pair and count this as the same hit

def main():
	"""Main function"""

	### Read in the resistant mutations and store in dictionaries
	# PI
	pi_dictionary_file = args.pi
	pi_in_fp = open(pi_dictionary_file, "r")
	pi_lines = pi_in_fp.readlines()
	pi_mutations = dict()
	for pi_line in pi_lines:
		pi_split_list = pi_line.split("\t")
		if len(pi_split_list) != 3: # Check to see if the line has 3 tab delimited fields, throw error if not
			print "ERROR: there was a problem with the formatting of the PI dictionary file - please check the format and try again. There must be no empty lines and each line must be tab separated with 3 columns for mutation, mutation class and comment."
			sys.exit()
		else:
			mutation, CLASS, comment = pi_split_list
			# print "PI ", mutation, " -> ", CLASS
			# pi_mutations[mutation] = ""
			pi_mutations[mutation] = CLASS

	# NRTI
	nrti_dictionary_file = args.nrti
	nrti_in_fp = open(nrti_dictionary_file, "r")
	nrti_lines = nrti_in_fp.readlines()
	nrti_mutations = dict()
	for nrti_line in nrti_lines:
		nrti_split_list = nrti_line.split("\t")
		if len(nrti_split_list) != 3: # Check to see if the line has 3 tab delimited fields, throw error if not
			print "ERROR: there was a problem with the formatting of the NRTI dictionary file - please check the format and try again. There must be no empty lines and each line must be tab separated with 3 columns for mutation, mutation class and comment."
			sys.exit()
		else:
			mutation, CLASS, comment = nrti_split_list
			# print "NRTI ", mutation, " -> ", CLASS
			# nrti_mutations[mutation] = ""
			nrti_mutations[mutation] = CLASS

	# NNRTI
	nnrti_dictionary_file = args.nnrti
	nnrti_in_fp = open(nnrti_dictionary_file, "r")
	nnrti_lines = nnrti_in_fp.readlines()
	nnrti_mutations = dict()
	for nnrti_line in nnrti_lines:
		nnrti_split_list = nnrti_line.split("\t")
		if len(nnrti_split_list) != 3: # Check to see if the line has 3 tab delimited fields, throw error if not
			print "ERROR: there was a problem with the formatting of the NNRTI dictionary file - please check the format and try again. There must be no empty lines and each line must be tab separated with 3 columns for mutation, mutation class and comment."
			sys.exit()
		else:
			mutation, CLASS, comment = nnrti_split_list
			# print "NNRTI ", mutation, " -> ", CLASS
			# nnrti_mutations[mutation] = ""
			nnrti_mutations[mutation] = CLASS

	### Run BLASTx on the fasta file
	# TODO: need to check for blastx executable and throw error if not installed.
	if args.fasta is not None and args.xml is not None:
		sys.exit("Both a FASTA and an BLASTx XML file were specified in the input parameters, it is only possible to use one of these input types")
	elif args.fasta is not None:
		print "Running BLASTx search on fasta input file..."
		# Set fasta_infile variable from the fasta input file specified on the command line
		fasta_infile = args.fasta
		# Set basefile name that is used for creating output files
		base_filename = os.path.splitext(fasta_infile)[0]
		blast_outfile = base_filename + "_blastout.xml" # Set the blast xml output filename
		os.system("blastx -query " + fasta_infile + " -db all -out " + blast_outfile + " -outfmt 5 -max_hsps 1")
		print "BLASTx search complete..."
	else:
		if args.xml is None:
			sys.exit("No FASTA or BLASTx XML file were specified, you must use one of these input types")
		else:
			print "Using supplied BLASTx XML file"
			blast_outfile = args.xml
			base_filename = os.path.splitext(blast_outfile)[0]

	### Open the BLASTx xml search output file - this will be read in later using BioPython
	result=open(blast_outfile,"r")

	### Open the output files for write
	individual_counts_out = open(base_filename + "_counts.txt",'w')
	linked_list_out = open(base_filename + "_linked_list.txt",'w')
	linked_counts_out = open(base_filename + "_linked_counts.txt",'w')

	linkedcount = args.linkedcount
	print "Linked mutations with a count of", linkedcount, "or more will be reported..."

	### Initialise all the defaultdict for storing counts of mutations
	pr_counts = defaultdict(int)
	pr_counts_accessory = defaultdict(int)
	pr_counts_dosage = defaultdict(int)
	pr_counts_other = defaultdict(int)
	pr_counts_unk= defaultdict(int)
	pr_counts_total = defaultdict(int)
	rt_counts = defaultdict(int)
	rt_counts_other = defaultdict(int)
	rt_counts_unk = defaultdict(int)
	rt_counts_total = defaultdict(int)
	linked_mutation_counts = defaultdict(int)

	mutation_count = defaultdict(int)

	### Initialise counts
	total_overall_mutation_count = 0
	pr_all_count_total = 0
	rt_all_count_total = 0

	record_number = 0

	# Store the evalue cutoff setting from input arguments
	evalue = args.evalue

	linked_mutations = defaultdict(list) # Declare the linked mutation dictionary

	### Start parsing the blastx records
	print "Parsing BLASTx records..."
	blast_records= NCBIXML.parse(result)
	for blast_record in blast_records: # Loop through each sequence input blast search
		
		query_name = blast_record.query
		record_number += 1
		# print "\nQuery name -> ", query_name, record_number
		# BLAST object - http://biopython.org/DIST/docs/tutorial/Tutorial.html#fig:blastrecord
		# BLAST tutorial - https://github.com/mbourgey/Concordia_Workshop_Biopython

		# print ('query_id -> ', blast_record.query_id)
		# print ('query -> ', blast_record.query)

		for alignment in blast_record.alignments:
			
			# print('hit sequence:', alignment.accession)
			for hsp in alignment.hsps:
					# TODO: Make sure that the sequence that's aligned isn't also hit in another hsp to avoid duplicates *** How to do this?? Possible use the -max_target_seqs setting in blastx to only keep one aligned sequence?
					if hsp.expect < evalue:
						# print('****Alignment****')
						# print('sequence:', alignment.title)
						# print('sequence:', alignment.accession) 
						# print('ID:', alignment.accession)
						# print('length:', alignment.length)
						# print('score:', hsp.score)
						# print('gaps:', hsp.gaps)
						# print('e value:', hsp.expect)
						# print(hsp.query)
						# print(hsp.match)
						# print(hsp.sbjct)
						# print('query start:',hsp.query_start)
						# print('query end:',hsp.query_start)
						# print('subject start:',hsp.sbjct_start)
						# print('subject end:',hsp.sbjct_end)
						# print('align length:',hsp.align_length)
						# print('score:',hsp.score)
						# print('expect:',hsp.expect)

						# region_title = str(alignment.title)
						region_title = str(alignment.accession)

						c = 0
						query = list(hsp.query)
						subject = list(hsp.sbjct)
						subject_position = hsp.sbjct_start
						for char in hsp.match: # Loop through each hsp (should only be the top one since this is set in the initial blastx search done in the script)
							# print c, char, "->", subject_position
							if char == " " or not char.isalpha(): # char == "+" # There's a mismatch between query and subject for this AA so it's a mutation
								mutation_string = str(subject[c]) + str(subject_position) + str(query[c]) # Create the mutation string
								total_overall_mutation_count += 1 # Increment the total count of mutations
								# print mutation_string
								if re.match('rt.*', region_title) is not None: # This hit is a match to the rt region
									rt_counts_total[subject_position] += 1
									rt_all_count_total += 1
									if mutation_string in nrti_mutations: # Check if this mutation is listed in the NRTI dictionary of mutations for this region
										# print "NRTI -> ", mutation_string, " ---> ", nrti_mutations[mutation_string], " ---> ", query_name
										mutation_count[query_name + "\t" + "NRTI_" + nrti_mutations[mutation_string]] += 1 # Increment the mutation count for this sequence hit to find linked mutations
										linked_mutations[query_name + "\t" + "NRTI_"  + nrti_mutations[mutation_string]].append(mutation_string) # Add the mutation to the defaultdict so that when looping through the mutation_count you can also print the mutations and find linked one
										# print "ADDING ", mutation_string
										rt_counts[str(subject_position) + "\t" + mutation_string] += 1 # Increment the count for this mutation to get the overall mutation count for the sample
									elif mutation_string in nnrti_mutations: # Check if this mutation is listed in the NNRTI dictionary of mutations for this region
										# print "NNRTI -> ", mutation_string
										mutation_count[query_name + "\t" + "NNRTI_" + nnrti_mutations[mutation_string]] += 1
										linked_mutations[query_name + "\t" + "NNRTI_" + nnrti_mutations[mutation_string]].append(mutation_string)
										rt_counts[str(subject_position) + "\t" + mutation_string] += 1
									else: # Mutation is not listed in either dictionary for RT region so it's an unknown mutation
										if args.unk is True:
											mutation_count[query_name + "\t" + "RT_UNK"] += 1
											linked_mutations[query_name + "\t" + "RT_UNK"].append(mutation_string)
											rt_counts_unk[str(subject_position) + "\t" + mutation_string] += 1

								if re.match('pr.*', region_title) is not None:
									pr_counts_total[subject_position] += 1
									pr_all_count_total += 1
									if mutation_string in pi_mutations: # Check if this mutation is listed in the dictionary of mutations for this region
										# print "PR -> ", mutation_string
										mutation_count[query_name + "\t" + "PR_" + pi_mutations[mutation_string]] += 1
										linked_mutations[query_name + "\t" + "PR_" + pi_mutations[mutation_string]].append(mutation_string)
										# print "NOTDEFINED PR ADDING -> " + mutation_string + " to ---> " + pi_mutations[mutation_string] + " ---> " + str(linked_mutations[query_name + "\t" + "PR_" + pi_mutations[mutation_string]])
										# if mutation_string == "":
											# print "NOTDEFINED - *********--> " + mutation_string
										pr_counts[str(subject_position) + "\t" + mutation_string] += 1
									else:
										if args.unk is True:
											# print "PR_UNK -> ", mutation_string
											mutation_count[query_name + "\t" + "PR_UNK"] += 1
											linked_mutations[query_name + "\t" + "PR_UNK"].append(mutation_string)
											pr_counts_unk[str(subject_position) + "\t" + mutation_string] += 1
							c += 1
							subject_position += 1

	# print "Finished parsing BLASTx records..."

	for region_key, region_count in mutation_count.iteritems(): # Iterate through the mutation count for this hit and also print the defaultdict value for the array of mutations
		if region_count >= linkedcount: # Only print where there's more than one mutation
			# print "count -> ", region_count, " region_key -> ", region_key
			# linked_mutations_list_sorted = linked_mutations[region_key]
			# linked_mutations_list_sorted.sort(key=natural_keys)
			# linked_mutations_list_sorted = sorted(linked_mutations[region_key], key=natural_key)
			# print "BEFORE NOTDEFINED -> " + str(linked_mutations[region_key])
			linked_mutations_list_sorted = natsort.natsorted(linked_mutations[region_key])
			# linked_mutations_list_sorted = sorted(linked_mutations[region_key])
			# linked_mutations_string = ",".join(linked_mutations_list_sorted)
			linked_mutations_string = ",".join(linked_mutations[region_key])
			# print "AFTER NOTDEFINED -> " + str(linked_mutations_string)
			# if str(linked_mutations_list_sorted) == "[]":
				# print "NOTDEFINED -> " + region_key + "\t" + str(region_count) + "\t -----> " + str(linked_mutations_list_sorted) + " ---> " + str(linked_mutations[region_key])
			# linked_list_out.write("LINKED" + "\t" + region_key + "\t" + str(region_count) + "\t" + str(linked_mutations_list_sorted) + " -> without sorting -> " +  linked_mutations_string + " ------> " + str(linked_mutations[region_key]) + "\n")

			linked_list_out.write(region_key + "\t" + str(region_count) + "\t" + str(linked_mutations_list_sorted) + "\n")

			# linked_list_out.write("LINKED" + "\t" + region_key + "\t" + str(region_count) + "\t" + str(linked_mutations_list_sorted) + " -> without sorting -> " +  linked_mutations_string + " ------> " + str(linked_mutations[region_key]) + "\n")
			# print "TESTLINKED" + "\t" + region_key + "\t" + str(region_count) + "\t -----> " + str(linked_mutations_list_sorted)

			tmp_split = region_key.split("\t")
			region_alone = tmp_split[-1]
			
			linked_mutation_counts[region_alone + "\t" + linked_mutations_string] += 1
			# print "region alone -> ", region_alone, " | linked_mutations_string -> ", linked_mutations_string, " --> ", linked_mutation_counts[region_alone + "\t" + linked_mutations_string]
				# print "LINKED_MUTATION_COUNT -> " + region_alone + "\t" + linked_mutations_string + " -> " + str(linked_mutation_counts[region_alone + "\t" + linked_mutations_string])

	# print "PR occurence counts:"
	for pr_key, pr_count in sorted(pr_counts.iteritems(), key=lambda (k,v): (v,k)):
		# if pr_count > 1:
		tmp_split = pr_key.split("\t") # split key on tabs
		position = tmp_split[0] # get the amino acid position
		total = pr_counts_total[int(position)] # find the total number of mutations for this amino acid position
		# percentage = (float(pr_count) / float(total) ) * 100 # calculate the percentage frequency of this mutation at this position
		percentage = (float(pr_count) / pr_all_count_total ) * 100 # calculate the percentage frequency of this mutation at this position
		# print pr_key + "\t" + str(pr_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count"
		# individual_counts_out.write(pr_key + "\t" + str(pr_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count" + "\n")
		individual_counts_out.write(pr_key + "\t" + str(pr_count) + "\t" + str(pr_all_count_total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count" + "\n")

	if args.unk is True:
		# print "PR other occurence counts:"
		for pr_key, pr_count in sorted(pr_counts_unk.iteritems(), key=lambda (k,v): (v,k)):
			# if pr_count_other > 1:
			tmp_split = pr_key.split("\t")
			position = tmp_split[0]
			total = pr_counts_total[int(position)]
			# percentage = (float(pr_count) / float(total) ) * 100
			percentage = (float(pr_count) / pr_all_count_total ) * 100
			# print pr_key + "\t" + str(pr_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count_other"
			# individual_counts_out.write(pr_key + "\t" + str(pr_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count_unk" + "\n")
			individual_counts_out.write(pr_key + "\t" + str(pr_count) + "\t" + str(pr_all_count_total) + "\t" + format(percentage, ".2f") + "\t" + "PR_count_unk" + "\n")

	# print "RT occurence counts:"
	for rt_key, rt_count in sorted(rt_counts.iteritems(), key=lambda (k,v): (v,k)):
		# if rt_count > 1:
		tmp_split = rt_key.split("\t")
		position = tmp_split[0]
		total = rt_counts_total[int(position)]
		# percentage = (float(rt_count) / float(total) ) * 100
		percentage = (float(rt_count) / rt_all_count_total ) * 100
		# print rt_key + "\t" + str(rt_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count"
		# individual_counts_out.write(rt_key + "\t" + str(rt_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count" + "\n")
		individual_counts_out.write(rt_key + "\t" + str(rt_count) + "\t" + str(rt_all_count_total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count" + "\n")

	if args.unk is True:
		# print "RT other occurence counts:"
		for rt_key, rt_count in sorted(rt_counts_unk.iteritems(), key=lambda (k,v): (v,k)):
			# if rt_count_other > 1:
			tmp_split = rt_key.split("\t")
			position = tmp_split[0]
			total = rt_counts_total[int(position)]
			# percentage = (float(rt_count) / float(total) ) * 100
			percentage = (float(rt_count) / rt_all_count_total ) * 100
			# print rt_key + "\t" + str(rt_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count_other"
			# individual_counts_out.write(rt_key + "\t" + str(rt_count) + "\t" + str(total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count_unk" + "\n")
			individual_counts_out.write(rt_key + "\t" + str(rt_count) + "\t" + str(rt_all_count_total) + "\t" + format(percentage, ".2f") + "\t" + "RT_count_unk" + "\n")


	# Loop through the linked mutation counts and output to final file
	for key, count in sorted(linked_mutation_counts.iteritems(), key=lambda (k,v): (v,k)):
		linked_counts_out.write(key + "\t" + str(count) + "\n")

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	'''
	alist.sort(key=natural_keys) sorts in human order
	http://nedbatchelder.com/blog/200712/human_sorting.html
	(See Toothy's implementation in the comments)
	'''
	return [ atoi(c) for c in re.split('(\d+)', text) ]

if __name__ == '__main__':
	sys.exit(main())
