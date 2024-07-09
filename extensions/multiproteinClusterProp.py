#!/usr/bin/env python3

import os
import glob
import pandas as pd
import argparse

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
	parser.add_argument('-i', '--bulk-clust-dir',  help='Directory created from bulk clustering script.', required=True)
	parser.add_argument('-o', '--output-dir', default="Network Proportions", help='Name of directory to output files with protein proportions of networks')
	
	args = parser.parse_args()

	if not os.path.exists(args.output_dir):
		os.mkdir(args.output_dir)
	else:
		print(f"Warning: The directory \"{args.output_dir}\" already exists. Files may be overwritten.")

	# loop through each species
	for spec_dir in glob.glob(os.path.join(args.bulk_clust_dir, '*')):
		if os.path.isdir(spec_dir) and len(glob.glob(os.path.join(spec_dir, 'multiprotein_networks'))) > 0:
			if not os.path.exists(os.path.join(args.output_dir, os.path.basename(spec_dir))):
				os.mkdir(os.path.join(args.output_dir, os.path.basename(spec_dir)))

			# loop through each network
			for network in glob.glob(os.path.join(spec_dir, 'multiprotein_networks', 'network*')):
				# get proportions from summary file
				sum_file = glob.glob(os.path.join(network, '*_summary_stats.tsv'))[0]
				summary_df = pd.read_csv(sum_file, sep="\t")
				clust_cols = [x for x in summary_df.columns if "_Cluster_Count" in x]

				clust_props = dict(zip(list(range(1, len(clust_cols) + 1)), [0]*len(clust_cols)))

				for idx, row in summary_df.iterrows():
					in_clust_count = sum([1 for col in clust_cols if row[col] != 0])
					clust_props[in_clust_count] += 1

				for clust_count in clust_props.keys():
					clust_props[clust_count] = round(clust_props[clust_count] / len(summary_df.index), 3)

				prop_df = pd.DataFrame(list(clust_props.items()), columns=["Protein Cluster Count", "Proportion"])

				prop_df.to_csv(os.path.join(args.output_dir, os.path.basename(spec_dir), f"{os.path.basename(network)}_protein_cluster_count_proportions"), sep="\t", index=False)


if __name__ == "__main__":
	main()