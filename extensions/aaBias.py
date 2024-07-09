#!/usr/bin/env python3

import os
import glob
import pandas as pd
import argparse
from collections import defaultdict

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
	parser.add_argument('-i', '--enriched-peptide-dir',  help='Directory with enriched petide files for input', required=True)
	parser.add_argument('-p', '--peptide-file',  help='Name of fasta file containing aa peptides of interest. These will generally be peptides that are contained in a particular assay.', required=True)
	parser.add_argument('-o', '--output-file', default="aa_bias.tsv", help='Name of .tsv to output file with AA bias data')
	
	args = parser.parse_args()

	pep_2_counts = defaultdict(lambda: defaultdict(int))
	pep_file_props = defaultdict(float)

	# get props for peptide file
	with open( args.peptide_file ) as pep_file:
		# get counts 
		for line in pep_file.readlines():
			line = line.strip()
			if line.startswith('>'):
				curr_pep_name = line[1:]
			else:
				for let in line:
					pep_file_props[let] += 1
					pep_2_counts[curr_pep_name][let] += 1

	# get proportions
	total_counts = sum(pep_file_props.values())
	for let, count in pep_file_props.items():
		pep_file_props[let] = count / total_counts

	# get props for each set
	out_data = list()
	for enr_file in glob.glob(os.path.join(args.enriched_peptide_dir, "*enriched.txt")):
		pep_set_props = defaultdict(float)
		peptide_set = list()
		aa_bias = 0

		with open(enr_file) as file:
			for line in file.readlines():
				line = line.strip()
				if line:
					peptide_set.append(line)

		# get counts
		for pep in peptide_set:
			for let in pep_2_counts[pep].keys():
				pep_set_props[let] += pep_2_counts[pep][let]

		# get proportions
		total_counts = sum(pep_set_props.values())
		for let, count in pep_set_props.items():
			pep_set_props[let] = count / total_counts

		for let in pep_set_props.keys():
			aa_bias += abs(pep_file_props[let] - pep_set_props[let])

		out_data.append( ( os.path.basename(enr_file), len(peptide_set), round(aa_bias, 3) ) )

	out_df = pd.DataFrame(out_data, columns=["Filename", "Peptide Count", "AA Bias"]).sort_values(by="AA Bias")
	out_df.to_csv(args.output_file, sep="\t", index=False)

if __name__ == "__main__":
	main()