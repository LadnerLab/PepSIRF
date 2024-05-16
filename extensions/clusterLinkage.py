#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(description="Generates linkage scores between each cluster.")
	
	# Required arguments
	parser.add_argument("-i", "--input-dir", type=str, metavar="", required=True, help="Directory containing cluster files. Summary files will be ignored.")
	parser.add_argument("-p", "--cluster-prefix", type=str, metavar="", required=True, help="Prefix of cluster files to use from the input directory."	
												"Example: If cluster file is \"clusters_Filoviridae_GP.fasta.tsv\", then the prefix is \"clusters_Filoviridae_\".")
	parser.add_argument("-m", "--metadata", type=str, metavar="", required=True, help="Distance thresholds to use for hierarchical clustering. Multiple values may be provided, all of which should be between 0 and 1.")
	parser.add_argument("-o", "--output-dir", type=str, metavar="", required=True, help="Directory to save cluster files, each file will be each distance threshold.")

	# Optional arguments
	parser.add_argument("--cluster-suffix", default=".fasta.tsv", type=str, metavar="", required=False, help="Suffix of cluster files to use from the input directory.")
	parser.add_argument("--seq-name-header", default="SequenceName", type=str, metavar="", required=False, help="Header for sequence name column in the metadata file.")
	parser.add_argument("--nt-accession-header", default="NCBIaccession-NT", type=str, metavar="", required=False, help="Header for nucleotide accession number column in the metadata file.")
	parser.add_argument("--nt-accession-delim", default=",", type=str, metavar="", required=False, help="Delimiter for multiple nucleotide accession numbers associated with a single sequence in the metadata file.")

	args=parser.parse_args()

	find_linkage_scores(
		cluster_dir = args.input_dir,
		cluster_prefix = args.cluster_prefix,
		metadata = args.metadata,
		cluster_suffix = args.cluster_suffix,
		seq_header = args.seq_name_header,
		nt_header = args.nt_accession_header,
		nt_delim = args.nt_accession_delim,
		output_dir = args.output_dir
		)

#------------------------------------------

def find_linkage_scores(
	cluster_dir: str,
	cluster_prefix: str,
	metadata: str,
	cluster_suffix: str,
	seq_header: str,
	nt_header: str,
	nt_delim: str,
	output_dir: str
	)->None:

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")
	
	# create metadata dataframe
	mdDf = pd.read_csv(metadata, sep="\t", index_col=seq_header)

	# create dictionary with distance threshold: cluster id: cluster number: {set of nucleotide accessions}
	NT_dict = create_NT_accession_dict(mdDf, nt_header, cluster_dir, cluster_prefix, cluster_suffix)

	# calulate linkage scores and create output dataframe for each dist_thresh
	for dist_thresh in NT_dict.keys():
		out_data = list()

		# loop through each id 1
		for id_1 in NT_dict[dist_thresh]:
			# loop throuch each cluster 1
			for clust_1 in NT_dict[dist_thresh][id_1]:
				# loop through each following id 2
				for id_2 in NT_dict[dist_thresh]:
					if id_1 < id_2:
						# loop through each cluster 2
						for clust_2 in NT_dict[dist_thresh][id_2]:
							# initialize linkage score
							linkage_score = 0

							# convert entire nucleotide accession numbers in cluster 2 into set 2
							NT_set_2 = set()
							for NT_accession in NT_dict[dist_thresh][id_2][clust_2]:
								NT_set_2.update( set(NT_accession.split(nt_delim)) )

							# loop through each nucleotide accession number in cluster 1
							for NT_accession in NT_dict[dist_thresh][id_1][clust_1]:
								# test if pair exists with set 2
								for curr_accession in NT_accession.split(nt_delim):
									if curr_accession in NT_set_2:
										# increment linkage score
										linkage_score += 1
										break

							# do it the other way

							# convert entire nucleotide accession numbers in cluster 1 into set 1
							NT_set_1 = set()
							for NT_accession in NT_dict[dist_thresh][id_1][clust_1]:
								NT_set_1.update( set(NT_accession.split(nt_delim)) )

							# loop through each nucleotide accession number in cluster 2
							for NT_accession in NT_dict[dist_thresh][id_2][clust_2]:
								# test if pair exists with set 1
								for curr_accession in NT_accession.split(nt_delim):
									if curr_accession in NT_set_1:
										# increment linkage score
										linkage_score += 1
										break

							# divide linkage score by 2 (get average)
							linkage_score /= 2

							# add data
							if linkage_score > 0:
								out_data.append( (id_1, clust_1, id_2, clust_2, linkage_score) )

		outDf = pd.DataFrame( out_data, columns=["p1", "c1", "p2", "c2", "linkageScore"] )

		outDf.to_csv(f"{output_dir}/{cluster_prefix}{dist_thresh}_linkage_scores.tsv", sep="\t", index=False)

def create_NT_accession_dict(mdDf, nt_header, cluster_dir, cluster_prefix, cluster_suffix):
	# collect cluster files (sort alphabetically)
	clust_files = sorted(glob.glob(f"{cluster_dir}/{cluster_prefix}*{cluster_suffix}"))

	# create SequenceName to NCBIaccession-NT dictionary
	seq_2_NT = mdDf[nt_header].to_dict()

	NT_dict = defaultdict(dict)

	for file in clust_files:
		# read file and remove sequences that are not associated with a cluster
		clustDf = pd.read_csv(file, sep="\t", index_col="Sequence").dropna()

		# extract id from filename
		file_id = file[len(f"{cluster_dir}/{cluster_prefix}"):len(file)-len(cluster_suffix)]

		# loop through each distance threshhold
		for dist_thresh in list(clustDf.columns):
			# create dictionary of set of unique nucleotide accession numbers assigned to each cluster
			clust_2_NTs = defaultdict(set)

			seq_2_clust = clustDf[dist_thresh].to_dict()
			for seq, clust in seq_2_clust.items():
				# check if sequence is assigned to a cluster
				clust_2_NTs[ int(clust) ].add(seq_2_NT[ seq ])

			NT_dict[dist_thresh][file_id] = clust_2_NTs

	return NT_dict

#------------------------------------------
if __name__ == "__main__":
	main();