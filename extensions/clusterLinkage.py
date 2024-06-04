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
NA = "NULL"

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
	parser.add_argument("--linkage-cols", nargs="+", default=["NCBIaccession-NT"], type=str, metavar="", required=False, help="Header of columns in the metadata file t consider for linkage. "
																											"Any given sequence only contributes one point to a linkage score")
	parser.add_argument("--col-val-delim", default=",", type=str, metavar="", required=False, help="Delimiter for multiple values of a column associated with a single sequence in the metadata file.")
	parser.add_argument("--make-network-vis", required=False, default=False, action="store_true", help="If provided, network visualization will be created for each threshold. The size of the nodes are based on the "
																										"number of sequences in the cluster and the size of the edges are based on the normalized linkage score.")
	parser.add_argument("--vis-output-dir", type=str, metavar="", required=False, default="networks_visualizations", help="Directory to save graphs if graphs are made.")
	parser.add_argument("--vis-seed", type=int, metavar="", required=False, default="1234", help="Seed used to generate network visualization. Changes node positions.")



	args=parser.parse_args()

	find_linkage_scores(
		cluster_dir = args.input_dir,
		cluster_prefix = args.cluster_prefix,
		metadata = args.metadata,
		cluster_suffix = args.cluster_suffix,
		seq_header = args.seq_name_header,
		linkage_cols = args.linkage_cols,
		col_val_delim = args.col_val_delim,
		make_net_vis = args.make_network_vis,
		vis_output_dir = args.vis_output_dir,
		vis_seed = args.vis_seed,
		output_dir = args.output_dir
		)

#------------------------------------------

def find_linkage_scores(
	cluster_dir: str,
	cluster_prefix: str,
	metadata: str,
	cluster_suffix: str,
	seq_header: str,
	linkage_cols: list,
	col_val_delim: str,
	make_net_vis: bool,
	vis_output_dir: str,
	vis_seed: int,
	output_dir: str
	)->None:

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")
	if make_net_vis:
		if not os.path.exists(vis_output_dir):
			os.mkdir(vis_output_dir)
		elif vis_output_dir != output_dir:
			print(f"Warning: The directory \"{vis_output_dir}\" already exists. Graphs may be overwritten.")
	
	# create metadata dataframe
	mdDf = pd.read_csv(metadata, sep="\t", index_col=seq_header)

	# create dictionary with distance threshold: cluster id: cluster number: sequence: [nucleotide accession, isolate, ...]
	data_dict = create_data_dict(mdDf, linkage_cols, cluster_dir, cluster_prefix, cluster_suffix)

	# create dictionary with distance threshold: cluster id: [{set of nucleotide accessions for the entire ID}, {set of all isolate}, ...]
	id_vals_dict = defaultdict(dict)
	for dist_thresh in data_dict:
		for id_key in data_dict[dist_thresh]:
			# create NT accession "cluster" for the entire id
			vals_list = list()
			for val_idx in range(len(linkage_cols)):
				value_set = set()
				for clust in data_dict[dist_thresh][id_key]:
					for seq in data_dict[dist_thresh][id_key][clust]:
						value_set.add(data_dict[dist_thresh][id_key][clust][seq][val_idx])
				vals_list.append(value_set)
			id_vals_dict[dist_thresh][id_key] = vals_list

	# calulate linkage scores and create output dataframe for each dist_thresh
	for dist_thresh in data_dict.keys():
		out_data = list()

		# loop through each id 1
		for id_1 in data_dict[dist_thresh]:
			# loop throuch each cluster 1
			for clust_1 in data_dict[dist_thresh][id_1]:
				# loop through each following id 2
				for id_2 in data_dict[dist_thresh]:
					if id_1 < id_2:
						# loop through each cluster 2
						for clust_2 in data_dict[dist_thresh][id_2]:
							linkage_score_1 = calc_linkage_score_from_c1_to_c2(
														clust_1=data_dict[dist_thresh][id_1][clust_1], 
														clust_2=data_dict[dist_thresh][id_2][clust_2], 
														col_val_delim=col_val_delim,
														num_metacols=len(linkage_cols)
														)

							if linkage_score_1 > 0:
								max_linkage_score_1 = calc_linkage_score_from_c1_to_c2(
															clust_1=data_dict[dist_thresh][id_1][clust_1], 
															clust_2=id_vals_dict[dist_thresh][id_2], 
															col_val_delim=col_val_delim,
															num_metacols=len(linkage_cols),
															findMax=True
														)
								norm_linkage_score = linkage_score_1 / max_linkage_score_1
								out_data.append( (id_1, clust_1, id_2, clust_2, norm_linkage_score) )

							# do it the otherway
							linkage_score_2 = calc_linkage_score_from_c1_to_c2(
														clust_1=data_dict[dist_thresh][id_2][clust_2], 
														clust_2=data_dict[dist_thresh][id_1][clust_1], 
														col_val_delim=col_val_delim,
														num_metacols=len(linkage_cols)
														)

							if linkage_score_2 > 0:
								max_linkage_score_2 = calc_linkage_score_from_c1_to_c2(
															clust_1=data_dict[dist_thresh][id_2][clust_2], 
															clust_2=id_vals_dict[dist_thresh][id_1], 
															col_val_delim=col_val_delim,
															num_metacols=len(linkage_cols),
															findMax=True
															)
								norm_linkage_score = linkage_score_2 / max_linkage_score_2
								out_data.append( (id_2, clust_2, id_1, clust_1, norm_linkage_score) )

		outDf = pd.DataFrame( out_data, columns=["p1", "c1", "p2", "c2", "normalizedLinkageScore"] )

		outDf.to_csv(f"{output_dir}/{cluster_prefix}{dist_thresh}_linkage_scores.tsv", sep="\t", index=False)

		create_network_visualization( dist_thresh, data_dict, outDf, cluster_prefix, vis_output_dir, vis_seed, output_dir, make_net_vis)


# calculate linkage score based on nucleotide accession numbers between clusters
def calc_linkage_score_between_clusts(clust_1, clust_2, col_val_delim):
	return (calc_linkage_score_from_c1_to_c2(clust_1, clust_2, col_val_delim) + 
		calc_linkage_score_from_c1_to_c2(clust_2, clust_1, col_val_delim)) / 2

# if findMax is True, clust 2 should be an array of sets of ALL values associated with an id
def calc_linkage_score_from_c1_to_c2(clust_1, clust_2, col_val_delim, num_metacols, findMax=False):
	# initialize linkage score
	linkage_score = 0
	# keep track of sequence contributed to linkage score
	contributed_sequences = set()

	# loop through each linkage columns
	for col_num in range( num_metacols ):
		val_set_2 = set()
		if findMax:
			for cur_val in clust_2[col_num]:
				val_set_2.update( set(cur_val.split(col_val_delim)) )
		else:
			for sequence in clust_2:
				val_set_2.update( set(clust_2[sequence][col_num].split(col_val_delim)) )

		# loop through each nucleotide accession number in cluster 1
		for sequence in clust_1:
			if sequence not in contributed_sequences:
				# test if pair exists with set 2
				for cur_val in clust_1[sequence][col_num].split(col_val_delim):
					if cur_val != NA and cur_val in val_set_2:
						# increment linkage score
						linkage_score += 1
						contributed_sequences.add(sequence)
						break

	return linkage_score

def create_data_dict(mdDf, linkage_cols, cluster_dir, cluster_prefix, cluster_suffix):
	# collect cluster files (sort alphabetically)
	clust_files = sorted(glob.glob(f"{cluster_dir}/{cluster_prefix}*{cluster_suffix}"))

	# create SequenceName to column values dictionary
	seq_2_values = mdDf[linkage_cols].fillna(NA).astype(str).transpose().to_dict('list')

	data_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

	for file in clust_files:
		# read file and remove sequences that are not associated with a cluster
		clustDf = pd.read_csv(file, sep="\t", index_col="Sequence").dropna()

		# extract id from filename
		file_id = file[len(f"{cluster_dir}/{cluster_prefix}"):len(file)-len(cluster_suffix)]

		# loop through each distance threshhold
		for dist_thresh in list(clustDf.columns):
			seq_2_clust = clustDf[dist_thresh].to_dict()

			for seq, clust in seq_2_clust.items():
				# check if sequence is assigned to a cluster, add values to dict
				data_dict[dist_thresh][file_id][int(clust)][seq] = seq_2_values[seq]		

	return data_dict

def create_network_visualization( dist_thresh, data_dict, outDf,  cluster_prefix, vis_output_dir, vis_seed, output_dir, make_net_vis):
	G = nx.MultiDiGraph()
	color_assigned = defaultdict()
	color_index = 0
	cluster_colors = list()

	# add nodes from first column clusters
	clusters = list(set((outDf["p1"] + "_" + outDf["c1"].astype(str)).to_list()))
	clusters.sort()
	cluster_sizes = list()
	for cluster in clusters:
		cluster = cluster.split("_")
		id_ = cluster[0]
		num = int(cluster[1])
		cluster_sizes.append(len(data_dict[dist_thresh][id_][num]) * 25)

		if id_ not in color_assigned.keys():
			color_assigned[id_] = COLOR_LIST[color_index]
			color_index += 1

		# assign color for id
		cluster_colors.append(color_assigned[id_])
	G.add_nodes_from(clusters)

	# add edges, thickness based on linkage scores
	for index, row in outDf.iterrows():
		c1 = f"{row['p1']}_{row['c1']}"
		c2 = f"{row['p2']}_{row['c2']}"
		G.add_edge(c1, c2, color="k", weight=row["normalizedLinkageScore"] * 2)

	if make_net_vis:
		'''
		comps = nx.strongly_connected_components(G)
		pos = {}
		for net_num, comp in enumerate(comps):
			subgraph = G.subgraph(comp)
			subgraph_pos = nx.spring_layout(subgraph)

			# Offset positions to space out clusters
			offset = np.array([net_num * 5, 0])
			for node in subgraph_pos:
				pos[node] = subgraph_pos[node] + offset
		'''
		weights = nx.get_edge_attributes(G,'weight').values()
		pos=nx.spring_layout(G, k=0.75, seed=vis_seed)
		fig, ax = plt.subplots(figsize = (20, 14))
		nx.draw_networkx_nodes(G, pos, ax=ax, node_size=cluster_sizes, node_color=cluster_colors)
		nx.draw_networkx_edges(G, pos, ax=ax, connectionstyle=f'arc3, rad = 0.25', arrows=True, width=list(weights))
		nx.draw_networkx_labels(G, pos, ax=ax)
		fig.savefig(f"{vis_output_dir}/{cluster_prefix}{dist_thresh}_visualization.png", bbox_inches='tight', dpi=300)

	#---------Summary Statistics-----------
	seq_data = list()

	orphaned_cluster_dict = defaultdict(set)
	for id_ in data_dict[dist_thresh].keys():
		for num in data_dict[dist_thresh][id_]:
			orphaned_cluster_dict[id_].add(num)

	counts_dict = defaultdict(list)
	net_stats = list()
	comps = nx.strongly_connected_components(G)
	for net_num, comp in enumerate(comps):
		# initialze dict for network num
		for id_ in data_dict[dist_thresh].keys():
			counts_dict[id_].append(0)

		# find average linkage score
		avg_linkage_score = 0
		edge_count = 0
		for i, clust_1 in enumerate(comp):
			clust_1_list = clust_1.split("_")
			id_1 = clust_1_list[0]
			num_1 = int(clust_1_list[1])

			orphaned_cluster_dict[id_1].remove(num_1)

			counts_dict[id_1][net_num] += 1

			connecting_edges = G.out_edges(clust_1)
			for j, clust_2 in enumerate(comp):
				if i != j:
					clust_2_list = clust_2.split("_")
					id_2 = clust_2_list[0]
					num_2 = int(clust_2_list[1])

					if( (clust_1, clust_2) in connecting_edges ):
						# print(f"{clust_1}, {clust_2}: {outDf[(outDf['p1']==id_1) & (outDf['c1']==num_1) & (outDf['p2']==id_2) & (outDf['c2']==num_2)]['normalizedLinkageScore'].squeeze()}")
						avg_linkage_score += outDf[(outDf['p1']==id_1) & (outDf['c1']==num_1) & (outDf['p2']==id_2) & (outDf['c2']==num_2)]['normalizedLinkageScore'].squeeze()
						edge_count += 1

		avg_linkage_score /= edge_count

		# find number of sequences
		seq_num = 0
		for i, clust in enumerate(comp):
			clust_list = clust.split("_")
			id_ = clust_list[0]
			num = int(clust_list[1])
			seq_num += len(data_dict[dist_thresh][id_][num])

			# assighn network number to sequence
			for seq in data_dict[dist_thresh][id_][num].keys():
				seq_data.append( (seq, net_num) )

		net_stats.append( (net_num, avg_linkage_score, seq_num) )
	
	# get orphaned clusters stats
	orphaned_data = defaultdict()
	# go through each orphaned cluster
	for orp_id in orphaned_cluster_dict.keys():
		for orp_num in orphaned_cluster_dict[orp_id]:
			seq_num = 0
			net_num += 1
			# add a row
			for id_ in data_dict[dist_thresh].keys():
				counts_dict[id_].append(0)
			# increment orphaned cluster
			counts_dict[orp_id][net_num] += 1

			for seq in data_dict[dist_thresh][orp_id][orp_num].keys():
				seq_data.append( (seq, net_num) )

			# app row to net stats
			net_stats.append( (net_num, 0, len(data_dict[dist_thresh][orp_id][orp_num])) )

	countsDf = pd.DataFrame.from_dict(counts_dict).rename_axis("MultiProteinCluster")
	statsDf = pd.DataFrame( net_stats, columns=["MultiProteinCluster", "AvgNormLinkageScore", "NumberOfSequences"]).set_index("MultiProteinCluster")
	statsDf = countsDf.merge(statsDf, on="MultiProteinCluster", how="right")
	for id_ in data_dict[dist_thresh].keys():
		statsDf.rename( columns = {id_:f"{id_}_Clusters"}, inplace = True)
	statsDf.to_csv(f"{output_dir}/{cluster_prefix}{dist_thresh}_summary_stats.tsv", sep="\t")

	seqDf = pd.DataFrame(seq_data, columns=["Sequence", "MultiProteinCluster"]).set_index("Sequence")
	seqDf.to_csv(f"{output_dir}/{cluster_prefix}{dist_thresh}_multiprotein_cluster_sequences.tsv", sep="\t")
	#--------------------------------------

#------------------------------------------
if __name__ == "__main__":
	main();