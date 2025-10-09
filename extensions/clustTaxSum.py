#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import glob
from collections import defaultdict
from natsort import natsorted

def main():
	parser = argparse.ArgumentParser(description="Generates taxomic summary from multiprotein cluster sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	# Required arguments
	parser.add_argument("-i", "--input", type=str, required=True, help="Directory containing a merge multiprotein cluster sequence file corresponding to each network. Can be generated using bulkClustering.py. A single filename can also be used.")
	parser.add_argument("-m", "--metadata", type=str, required=True, help="Tab-delimited metadata file which contains taxomic info")
	parser.add_argument("-o", "--output", type=str, required=True, help="Directory to save taxomic summary of each network. This directory will be created, if it doesn't already exist. If a single file is inputted, the output should be a single filename")

	# Optional arguments
	parser.add_argument("--tax-cols", nargs="+", default=["Species", "Genus"], type=str, required=False, help="Header of columns in the metadata file to generate summary data for.")

	args=parser.parse_args()

	gen_sum(
		input_dir = args.input,
		metadata = args.metadata,
		tax_cols = args.tax_cols,
		output_dir = args.output,
		)

###---------------End of main()-----------------------------------


def gen_sum(
	input_dir: str,
	metadata: str,
	tax_cols: list,
	output_dir: str
	) -> None:
	
	if not os.path.isdir(input_dir):
		files = [input_dir]

	else:
		files = glob.glob(os.path.join(input_dir, 'network_*_merged_multiprotein_cluster_sequences.tsv'))
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
		else:
			print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")

	# get columns from metadata
	metaDf = pd.read_csv(metadata, sep="\t").set_index("SequenceName")[tax_cols]

	# look through each input file
	for file in files:
		net_num = os.path.basename(file)[len("network_"):len(os.path.basename(file)) - len("_merged_multiprotein_cluster_sequences.tsv")]

		inDf = pd.read_csv(file, sep="\t")
		seq_clust_pairs = set(zip(inDf["Sequence"], inDf["MultiProteinCluster"]))

		# create dictionary of cluster and sequences
		clust_seqs_dict = defaultdict(set)
		for pair in seq_clust_pairs:
			clust_seqs_dict[pair[1]].add(pair[0])

		# get data for output
		outDf = pd.DataFrame([(x, len(clust_seqs_dict[x])) for x in natsorted(clust_seqs_dict.keys(), key=lambda key: key.lower())], columns=["MultiProteinCluster", "Sequence_Count"])
		for col in tax_cols:
			col_out_data = list()
			for clust in clust_seqs_dict.keys():
				col_counts = defaultdict()
				total_count = 0

				for seq in clust_seqs_dict[clust]:
					# add count to corresponding column
					col_name = metaDf.loc[[seq], col].item()

					if col_name not in col_counts.keys():
						col_counts[col_name] = 1
					else:
						col_counts[col_name] += 1

					total_count += 1

				# get name list (sorted by count)
				col_names = [k for k,v in sorted(col_counts.items(), key=lambda item: item[1], reverse=True)]
				# get total count
				total_name_count = len(col_names)
				# get proportion list
				prop_list = list()
				for name in col_names:
					prop_list.append(round(col_counts[name] / total_count, 3))

				col_out_data.append((clust, total_name_count, ",".join([str(x) for x in col_names]), ",".join([str(x) for x in prop_list])))

			colDf = pd.DataFrame( col_out_data, columns=["MultiProteinCluster", f"{col}_Count", f"{col}_List", f"{col}_Prop"] )
			outDf = outDf.merge(colDf, on="MultiProteinCluster", how="left")

		if not os.path.isdir(input_dir):
			outDf.to_csv(os.path.join(output_dir), sep="\t", index=False)
		else:
			outDf.to_csv(os.path.join(output_dir, f"network_{net_num}_taxomic_summary.tsv"), sep="\t", index=False)
		

if __name__ == "__main__":
	main()
