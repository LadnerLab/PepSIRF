#!/usr/bin/env python3

import os
import glob
import pandas as pd
import argparse

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
	parser.add_argument('-i', '--bulk-clust-dir', help='Directory created from bulk clustering script.', required=True)
	parser.add_argument('--prop-thresh', type=float, default=0.75, help='Proportion theshold to alert user if greater than this proportion of clusters contains less than or equal to a specified count of proteins', required=False)
	parser.add_argument('--count-thresh', type=int, default=1, help='Count threshold to alert user if less than or equal to this count of different proteins is contained in greater than a specified proportion of the clusters', required=False)
	parser.add_argument('-o', '--output-dir', default="Network Proportions", help='Name of directory to output files with protein proportions of networks')
	
	args = parser.parse_args()

	assert args.prop_thresh >= 0 and args.prop_thresh <= 1, \
		f"Warning: Proportion threshold is not a proportion"

	if not os.path.exists(args.output_dir):
		os.mkdir(args.output_dir)
	else:
		print(f"Warning: The directory \"{args.output_dir}\" already exists. Files may be overwritten.\n")

	# loop through each species
	for spec_dir in glob.glob(os.path.join(args.bulk_clust_dir, '*')):
		if os.path.isdir(spec_dir) and len(glob.glob(os.path.join(spec_dir, 'multiprotein_networks'))) > 0:

			# check for greater than 1 protein cluster:
			if len(glob.glob(os.path.join(spec_dir, '*.fasta'))) > 1:
				if not os.path.exists(os.path.join(args.output_dir, os.path.basename(spec_dir))):
					os.mkdir(os.path.join(args.output_dir, os.path.basename(spec_dir)))

				# loop through each network
				for network in glob.glob(os.path.join(spec_dir, 'multiprotein_networks', 'network*')):
					# get proportions from summary file
					sum_file = glob.glob(os.path.join(network, '*_summary_stats.tsv'))[0]
					summary_df = pd.read_csv(sum_file, sep="\t")
					clust_cols = [x for x in summary_df.columns if "_Cluster_Count" in x]
					out_data = list()

					clust_counts = dict(zip(list(range(1, len(clust_cols) + 1)), [0]*len(clust_cols)))

					for idx, row in summary_df.iterrows():
						in_clust_count = sum([1 for col in clust_cols if row[col] != 0])
						clust_counts[in_clust_count] += 1

					clust_thresh_prop = 0
					for clust_count in clust_counts.keys():
						prop = round(clust_counts[clust_count] / len(summary_df.index), 3)
						out_data.append((clust_count, prop, clust_counts[clust_count]))

						if clust_count <= args.count_thresh:
							clust_thresh_prop += prop

					if clust_thresh_prop > args.prop_thresh:
						print(f"In {os.path.basename(spec_dir)} {os.path.basename(network)}, over {args.prop_thresh * 100}% of clusters contain less than or equal to {args.count_thresh} protein(s)")

					prop_df = pd.DataFrame(out_data, columns=["Protein Cluster Count", "Proportion", "Count"])

					prop_df.to_csv(os.path.join(args.output_dir, os.path.basename(spec_dir), f"{os.path.basename(network)}_protein_cluster_count_proportions.tsv"), sep="\t", index=False)
			else:
				print(f"{os.path.basename(spec_dir)} skipped because it only contains 1 protein.")


if __name__ == "__main__":
	main()