#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

COLOR_LIST = ['#ff0000', '#deff0a', '#0aefff', '#be0aff', '#ff8700', '#a1ff0a', '#147df5', '#ffd300', '#0aff99', '#580aff']

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
	parser.add_argument("--make-network-vis", required=False, default=False, action="store_true", help="If provided, network visualization will be created for each threshold. The size of the nodes are based on the "
																										"number of sequences in the cluster and the size of the edges are based on the normalized linkage score.")
	parser.add_argument("--vis-output-dir", type=str, metavar="", required=False, default="networks_visualizations", help="Directory to save graphs if graphs are made.")

	args=parser.parse_args()

	find_linkage_scores(
		cluster_dir = args.input_dir,
		cluster_prefix = args.cluster_prefix,
		metadata = args.metadata,
		cluster_suffix = args.cluster_suffix,
		seq_header = args.seq_name_header,
		nt_header = args.nt_accession_header,
		nt_delim = args.nt_accession_delim,
		make_net_vis = args.make_network_vis,
		vis_output_dir = args.vis_output_dir,
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
	make_net_vis: bool,
	vis_output_dir: str,
	output_dir: str
	)->None:

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")
	if make_net_vis:
		if not os.path.exists(vis_output_dir):
			os.mkdir(vis_output_dir)
		else:
			print(f"Warning: The directory \"{vis_output_dir}\" already exists. Graphs may be overwritten.")
	
	# create metadata dataframe
	mdDf = pd.read_csv(metadata, sep="\t", index_col=seq_header)

	# create dictionary with distance threshold: cluster id: cluster number: {set of nucleotide accessions}
	nt_dict = create_NT_accession_dict(mdDf, nt_header, cluster_dir, cluster_prefix, cluster_suffix)
	
	# create dictionary with distance threshold: cluster id: {set of nucleotide accessions for the entire ID}
	id_nt_dict = defaultdict(dict)
	for dist_thresh in nt_dict:
		for id_key in nt_dict[dist_thresh]:
			# create NT accession "cluster" for the entire id
			nt_set = set()
			for clust in nt_dict[dist_thresh][id_key]:
				nt_set.update(nt_dict[dist_thresh][id_key][clust])
			id_nt_dict[dist_thresh][id_key] = nt_set

	# calulate linkage scores and create output dataframe for each dist_thresh
	for dist_thresh in nt_dict.keys():
		out_data = list()

		# loop through each id 1
		for id_1 in nt_dict[dist_thresh]:
			# loop throuch each cluster 1
			for clust_1 in nt_dict[dist_thresh][id_1]:
				# loop through each following id 2
				for id_2 in nt_dict[dist_thresh]:
					if id_1 < id_2:
						# loop through each cluster 2
						for clust_2 in nt_dict[dist_thresh][id_2]:
							linkage_score_1 = calc_linkage_score_from_c1_to_c2(
														nt_clust_1=nt_dict[dist_thresh][id_1][clust_1], 
														nt_clust_2=nt_dict[dist_thresh][id_2][clust_2], 
														nt_delim=nt_delim
														)

							if linkage_score_1 > 0:
								max_linkage_score_1 = calc_linkage_score_from_c1_to_c2(
															nt_clust_1=nt_dict[dist_thresh][id_1][clust_1], 
															nt_clust_2=id_nt_dict[dist_thresh][id_2], 
															nt_delim=nt_delim
														)
								norm_linkage_score = linkage_score_1 / max_linkage_score_1
								out_data.append( (id_1, clust_1, id_2, clust_2, norm_linkage_score) )

							linkage_score_1 = calc_linkage_score_from_c1_to_c2(
														nt_clust_1=nt_dict[dist_thresh][id_1][clust_1], 
														nt_clust_2=nt_dict[dist_thresh][id_2][clust_2], 
														nt_delim=nt_delim
														)

							if linkage_score_1 > 0:
								max_linkage_score_1 = calc_linkage_score_from_c1_to_c2(
															nt_clust_1=nt_dict[dist_thresh][id_1][clust_1], 
															nt_clust_2=id_nt_dict[dist_thresh][id_2], 
															nt_delim=nt_delim
														)
								norm_linkage_score = linkage_score_1 / max_linkage_score_1
								out_data.append( (id_1, clust_1, id_2, clust_2, norm_linkage_score) )

							# do it the otherway
							linkage_score_2 = calc_linkage_score_from_c1_to_c2(
														nt_clust_1=nt_dict[dist_thresh][id_2][clust_2], 
														nt_clust_2=nt_dict[dist_thresh][id_1][clust_1], 
														nt_delim=nt_delim
														)

							if linkage_score_2 > 0:
								max_linkage_score_2 = calc_linkage_score_from_c1_to_c2(
															nt_clust_1=nt_dict[dist_thresh][id_2][clust_2], 
															nt_clust_2=id_nt_dict[dist_thresh][id_1], 
															nt_delim=nt_delim
															)
								norm_linkage_score = linkage_score_2 / max_linkage_score_2
								out_data.append( (id_2, clust_2, id_1, clust_1, norm_linkage_score) )

		outDf = pd.DataFrame( out_data, columns=["p1", "c1", "p2", "c2", "normalizedLinkageScore"] )

		outDf.to_csv(f"{output_dir}/{cluster_prefix}{dist_thresh}_linkage_scores.tsv", sep="\t", index=False)

		if make_net_vis:
			create_network_visualization( dist_thresh, nt_dict, outDf, cluster_prefix, vis_output_dir)


# calculate linkage score based on nucleotide accession numbers between clusters
def calc_linkage_score_between_clusts(nt_clust_1, nt_clust_2, nt_delim):
	return (calc_linkage_score_from_c1_to_c2(nt_clust_1, nt_clust_2, nt_delim) + 
		calc_linkage_score_from_c1_to_c2(nt_clust_2, nt_clust_1, nt_delim)) / 2

def calc_linkage_score_from_c1_to_c2(nt_clust_1, nt_clust_2, nt_delim):
	# initialize linkage score
	linkage_score = 0

	# convert entire nucleotide accession numbers in cluster 2 into set 2
	NT_set_2 = set()
	for NT_accession in nt_clust_2:
		NT_set_2.update( set(NT_accession.split(nt_delim)) )

	# loop through each nucleotide accession number in cluster 1
	for NT_accession in nt_clust_1:
		# test if pair exists with set 2
		for curr_accession in NT_accession.split(nt_delim):
			if curr_accession in NT_set_2:
				# increment linkage score
				linkage_score += 1
				break

	return linkage_score

def create_NT_accession_dict(mdDf, nt_header, cluster_dir, cluster_prefix, cluster_suffix):
	# collect cluster files (sort alphabetically)
	clust_files = sorted(glob.glob(f"{cluster_dir}/{cluster_prefix}*{cluster_suffix}"))

	# create SequenceName to NCBIaccession-NT dictionary
	seq_2_NT = mdDf[nt_header].to_dict()

	nt_dict = defaultdict(dict)

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

			nt_dict[dist_thresh][file_id] = clust_2_NTs

	return nt_dict

def create_network_visualization( dist_thresh, nt_dict, outDf,  cluster_prefix, vis_output_dir ):
	G = nx.DiGraph()
	color_assigned = defaultdict()
	color_index = 0
	cluster_colors = list()

	# add nodes from first column clusters
	clusters = list(set((outDf["p1"] + "_" + outDf["c1"].astype(str)).to_list()))
	cluster_sizes = list()
	for cluster in clusters:
		cluster = cluster.split("_")
		id_ = cluster[0]
		num = int(cluster[1])
		cluster_sizes.append(len(nt_dict[dist_thresh][id_][num]) * 25)

		if id_ not in color_assigned.keys():
			color_assigned[id_] = COLOR_LIST[color_index]
			color_index += 1

		# assign color for id
		cluster_colors.append(color_assigned[id_])
	G.add_nodes_from(clusters)

	# add edges, thickness based on linkage scores
	weights = list()
	for index, row in outDf.iterrows():
		c1 = f"{row['p1']}_{row['c1']}"
		c2 = f"{row['p2']}_{row['c2']}"
		weights.append(row["normalizedLinkageScore"] * 2)
		G.add_edge(c1, c2)

	pos=nx.spring_layout(G, k=0.75, seed=5)
	fig, ax = plt.subplots(figsize = (20, 14))
	nx.draw_networkx_nodes(G, pos, ax=ax, node_size=cluster_sizes, node_color=cluster_colors)
	nx.draw_networkx_edges(G, pos, ax=ax, connectionstyle=f'arc3, rad = 0.25', arrows=True, width=weights)
	nx.draw_networkx_labels(G, pos, ax=ax)
	fig.savefig(f"{vis_output_dir}/{cluster_prefix}{dist_thresh}_visualization.png", bbox_inches='tight', dpi=300)

#------------------------------------------
if __name__ == "__main__":
	main();