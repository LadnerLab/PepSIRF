#!/usr/bin/env python

import argparse, os
import kmertools as kt
import fastatools as ft
import inout as io
from collections import defaultdict

# This script works like the link module in PepSIRF, but implements patterned kmers

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('--protein_file',  help='Name of fasta file containing target protein sequences.', required=True)
	reqArgs.add_argument('--peptide_file',  help='Name of fasta file containing aa peptides of interest. These will generally be peptides that are contained in a particular assay.', required=True)
	reqArgs.add_argument('--meta',  help='Taxonomic information for each protein contained in "--protein_file". Three comma-separated strings should be provided: 1) name of tab-delimited metadata file, 2) header for column containing protein sequence name and 3) header for column containing ID to be used in creating the linkage map.', required=True)

	parser.add_argument('-k', '--kmer_size', default=8, type=int, help='Number of positions at which there must be a match.')
	parser.add_argument('-s', '--span', default=10, type=int, help='Number of positions across which the match can span.')
	parser.add_argument('-o', '--output', default="link_output.tsv", help='Name for output linkage map.')
	
	args = parser.parse_args()

	# Generate all patterns and report the number being tested
	allPatts = kt.generatePatterns(args.kmer_size, args.span)
	print(f"Considering {len(allPatts)} patterns for generating linkage map.\n")
	
	# Read in metadata
	metaF, seqN, metaV = args.meta.split(",")
	metaD = io.fileDictHeaderLists(metaF, metaV, seqN)
	
	# Read in target sequences
	prot_fD = ft.read_fasta_dict(args.protein_file)
	# Read in peptide sequences
	pep_fD = ft.read_fasta_dict(args.peptide_file)
	
	outScoresD = defaultdict(list)
	
	# Determine total number of categories to process
	total = len(metaD)
	numProc=0
	
	# Step through each metadata category
	for each, nameL in metaD.items():
		print(each, len(nameL))
		# Step through each pattern of interest
		pepScoreD = defaultdict(int)
		for patt in allPatts:
			print(patt)
			targetKmers = {}
			# Step through each sequence for this group and add it's kmers to the dictionary
			for name in nameL:
				if name in prot_fD:
					for every in kt.patternedKmers(prot_fD[name], patt, outType="set", filter=["X", "B", "J", "Z"]):
						targetKmers[every] = ""
			
			this_pepScoreD = {}
			# Step through each peptide
			for pepN, pepS in pep_fD.items():
				pepKmers = kt.patternedKmers(pepS, patt, outType="set", filter=["X", "B", "J", "Z"])
				thisScore = sum([1 for kmer in pepKmers if kmer in targetKmers])
				this_pepScoreD[pepN] = thisScore
			#Update scores, if they're higher than for previous patterns
			for pepN, pepScore in this_pepScoreD.items():
				if pepScore > pepScoreD[pepN]:
					pepScoreD[pepN] = pepScore
		#for peptides with scores>0, add info to outScoresD
		for pepN, pepScore in pepScoreD.items():
			if pepScore > 0:
				outScoresD[pepN].append(f"{each}:{pepScore}")
		
		# Increment counter
		numProc+=1
		if numProc%1==0:
			print(numProc/total*100)
		
	#Write scores to output file
	with open(args.output, "w") as fout:
		fout.write("Peptide Name\tLinked Species IDs with counts\n")
		for pepN in pep_fD:
			if pepN in outScoresD:
				outStr = ",".join(outScoresD[pepN])
				fout.write(f"{pepN}\t{outStr}\n")


#----------------------End of main()



###------------------------------------->>>>	

if __name__ == "__main__":
	main()

