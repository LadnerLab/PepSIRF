#!/usr/bin/env python3

import pandas as pd
import argparse
import fastatools as ft
import numpy as np

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
	parser.add_argument('-i', '--fasta-file',  help='Directory with enriched petide files for input', required=True)
	parser.add_argument('-o', '--output-file', default="e_k_bias_out.tsv", help='Name of .tsv to output file with AA bias data')
	
	args = parser.parse_args()

	# get proportion of e's and k's for each peptide
	e_k_props = get_e_k_props(args.fasta_file)

	# get percentiles
	percentile_dict = get_percentiles(e_k_props)

	# create output df
	out_data = [(name, round(prop, 3), round(percentile_dict[name], 2)) for name, prop in e_k_props.items()]

	pd.DataFrame(out_data, columns=["CodeName", "e_k_Prop", "e_k_Percentile"]).to_csv(args.output_file, index=False, sep='\t')


# get proportion of e's and k's for each peptide
def get_e_k_props(fasta_file)->dict:
	e_k_props = dict()

	# get props for peptide file
	fasta_dict = ft.read_fasta_dict(fasta_file)

	# iterate through each sequence
	for name, seq in fasta_dict.items():
		e_k_count = 0

		# loop through each AA, get count of e and k
		for aa in seq:
			if aa.lower() == 'e' or aa.lower() =='k':
				e_k_count += 1

		# add proportion to dict
		e_k_props[name] = (e_k_count) / len(seq)

	return e_k_props

# get percentile of each peptide using its e and k proportion
def get_percentiles(e_k_props)->dict:
	# Calculate percentile of each peptide
	names = list(e_k_props.keys())
	all_props = np.array(list(e_k_props.values()))

	# Get unique values and array for mapping each original value to its corresponding index in the unique array
	unique_props, inverse_indices = np.unique(all_props, return_inverse=True)

	# Calculate percentiles based on the unique props
	percentile_ranks = np.linspace(0, 100, len(unique_props))

	# Map sorted values back to the original corresponding name
	percentile_dict = {names[i]: percentile_ranks[inverse_indices[i]] for i in range(len(names))}

	return percentile_dict


if __name__ == "__main__":
	main()