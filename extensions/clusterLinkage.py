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
	parser = argparse.ArgumentParser(description="Generates linkage scores between each cluster.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	# Required arguments
	parser.add_argument("-i", "--input-cluster-manifest", type=str, metavar="", required=True, help="Filepath of the protein clusters manifest file. The manifest file is tab delimited containing a "
                          											"header with the column names 'ProtID' and 'ClustersFile'. Each line specifies the protein and the filepath to the associated alignment file.")
	parser.add_argument("-m", "--metadata", type=str, metavar="", required=True, help="Distance thresholds to use for hierarchical clustering. Multiple values may be provided, all of which should be between 0 and 1.")
	parser.add_argument("-o", "--output-dir", type=str, metavar="", required=True, help="Directory to save output files.")

	# Optional arguments
	parser.add_argument("-t", "--thresh-matrix", type=str, metavar="", required=False, default="", help="Filepath to tab delimited file of the threshold matrix. Each column should represent the same protein id and the transposed rows "
																							"of the manifest file. Each row should be each seperate linkage network to create. The entries are the thresholds to use for each "
																							"protein id for each network to create. Optionally, an entry can be left blank and it will not be contributed to that network")
	parser.add_argument("--seq-name-header", default="SequenceName", type=str, metavar="", required=False, help="Header for sequence name column in the metadata file.")
	parser.add_argument("--linkage-cols", nargs="+", default=["NCBIaccession-NT"], type=str, metavar="", required=False, help="Header of columns in the metadata file t consider for linkage. "
																											"Any given sequence only contributes one point to a linkage score")
	parser.add_argument("--col-val-delim", default=",", type=str, metavar="", required=False, help="Delimiter for multiple values of a column associated with a single sequence in the metadata file.")
	parser.add_argument("--make-network-vis", required=False, default=False, action="store_true", help="If provided, network visualization will be created for each threshold. The size of the nodes are based on the "
																										"number of sequences in the cluster and the size of the edges are based on the normalized linkage score.")
	parser.add_argument("--vis-seed", type=int, metavar="", required=False, default="1234", help="Seed used to generate network visualization. Changes node positions.")



	args=parser.parse_args()

	find_linkage_scores(
		manifest_file = args.input_cluster_manifest,
		metadata = args.metadata,
		thresh_matrix = args.thresh_matrix,
		seq_header = args.seq_name_header,
		linkage_cols = args.linkage_cols,
		col_val_delim = args.col_val_delim,
		make_net_vis = args.make_network_vis,
		vis_seed = args.vis_seed,
		output_dir = args.output_dir
		)

#------------------------------------------

def find_linkage_scores(
	manifest_file: str,
	metadata: str,
	thresh_matrix,
	seq_header: str,
	linkage_cols: list,
	col_val_delim: str,
	make_net_vis: bool,
	vis_seed: int,
	output_dir: str
	)->None:
	# TODO: 2. let user specify threshold combination, if not use all together	

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")
	
	# create metadata dataframe
	mdDf = pd.read_csv(metadata, sep="\t", index_col=seq_header)

	# create dictionary with distance threshold: cluster id: cluster number: sequence: [nucleotide accession, isolate, ...]
	# also create thresh dictionary with: cluster id: [list of distance threshold to use for each index (network)]
	data_dict, thresh_dict = create_data_dict(mdDf, linkage_cols, manifest_file, thresh_matrix)

	# create dictionary with distance threshold: cluster id: sequence: [nucleotide accession, isolate, ...]
	id_vals_dict = defaultdict(lambda: defaultdict(dict))
	for dist_thresh in data_dict:
		for id_key in data_dict[dist_thresh]:
			for clust_num in data_dict[dist_thresh][id_key]:
				for seq in data_dict[dist_thresh][id_key][clust_num]:
					id_vals_dict[dist_thresh][id_key][seq] = data_dict[dist_thresh][id_key][clust_num][seq]

	# calulate linkage scores and create output dataframe for each dist_thresh
	for dist_idx in range(len(list(thresh_dict.values())[0])):
		net_out_dir = f"{output_dir}/network_{dist_idx +1}"
		if not os.path.exists(net_out_dir):
			os.mkdir(net_out_dir)

		out_data = list()

		ids = [prot_id for prot_id in thresh_dict.keys() if thresh_dict[prot_id][dist_idx] != NA]

		# loop through each id 1
		for id_1 in ids:
			thresh_1 = thresh_dict[id_1][dist_idx]
			# loop throuch each cluster 1
			for clust_1 in data_dict[thresh_1][id_1]:
				# loop through each following id 2
				for id_2 in ids:
					thresh_2 = thresh_dict[id_2][dist_idx]
					# loop through each cluster 2
					for clust_2 in data_dict[thresh_2][id_2]:
						if id_1 != id_2 or clust_1 != clust_2:
							linkage_score = calc_linkage_score_from_c1_to_c2(
														clust_1=data_dict[thresh_1][id_1][clust_1], 
														clust_2=data_dict[thresh_2][id_2][clust_2], 
														col_val_delim=col_val_delim,
														num_metacols=len(linkage_cols)
														)

							if linkage_score > 0:
								# remove sequences from max dict if they exist in clust 1
								clust_1_seqs = data_dict[thresh_1][id_1][clust_1].keys()
								max_id_2_dict = id_vals_dict[thresh_2][id_2].copy()
								for seq in clust_1_seqs:
									removed = max_id_2_dict.pop(seq, "Not Found")

								max_linkage_score = calc_linkage_score_from_c1_to_c2(
															clust_1=data_dict[thresh_1][id_1][clust_1], 
															clust_2=max_id_2_dict, 
															col_val_delim=col_val_delim,
															num_metacols=len(linkage_cols)
														)
								norm_linkage_score = linkage_score / max_linkage_score
								out_data.append( (id_1, clust_1, id_2, clust_2, norm_linkage_score) )

		outDf = pd.DataFrame( out_data, columns=["p1", "c1", "p2", "c2", "normalizedLinkageScore"] )

		outDf.to_csv(f"{net_out_dir}/network_{dist_idx +1}_linkage_scores.tsv", sep="\t", index=False)

		create_network( dist_idx, thresh_dict, ids, data_dict, outDf, vis_seed, net_out_dir, make_net_vis)

		# create a key linking network # to thresholds used
		thresh_key = []
		for prot_id in ids:
			thresh_key.append(thresh_dict[prot_id][dist_idx])
		thresh_key_df = pd.DataFrame( [thresh_key], columns=ids )
		thresh_key_df.to_csv(f"{net_out_dir}/network_{dist_idx +1}_thresh_key.tsv", sep="\t", index=False)


# calculate linkage score based on nucleotide accession numbers between clusters
def calc_linkage_score_between_clusts(clust_1, clust_2, col_val_delim):
	return (calc_linkage_score_from_c1_to_c2(clust_1, clust_2, col_val_delim) + 
		calc_linkage_score_from_c1_to_c2(clust_2, clust_1, col_val_delim)) / 2


# if findMax is True, clust 2 should be an array of sets of ALL values associated with an id
def calc_linkage_score_from_c1_to_c2(clust_1, clust_2, col_val_delim, num_metacols):
	# initialize linkage score
	linkage_score = 0
	# keep track of sequence contributed to linkage score
	contributed_sequences = set()

	# loop through each linkage columns
	for col_num in range( num_metacols ):
		val_set_2 = set()
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

def create_data_dict(mdDf, linkage_cols, manifest_file, thresh_matrix):
	# collect cluster files (sort alphabetically)
	id_2_file = pd.read_csv(manifest_file, sep="\t", index_col="ProtID").to_dict('index')
	ids = list(id_2_file.keys())

	thresh_dict = defaultdict(list)
	if thresh_matrix != "":
		thresh_data = pd.read_table(thresh_matrix, header=None).fillna(NA).values.tolist()

		for net in thresh_data:
			for i, prot_id in enumerate(ids):
				thresh_dict[prot_id].append(str(net[i]))

	# create SequenceName to column values dictionary
	seq_2_values = mdDf[linkage_cols].fillna(NA).astype(str).transpose().to_dict('list')

	data_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

	for prot_id in ids:
		# read file and remove sequences that are not associated with a cluster
		clustDf = pd.read_csv(id_2_file[prot_id]["ClustersFile"], sep="\t", index_col="Sequence").dropna()

		if thresh_matrix != "":
			dist_list = set(thresh_dict[prot_id])
		else:
			dist_list = set(clustDf.columns)
			thresh_dict[prot_id] = list(clustDf.columns)

		# loop through each distance threshold
		for dist_thresh in dist_list:
			if dist_thresh != NA:
				seq_2_clust = clustDf[dist_thresh].to_dict()

				for seq, clust in seq_2_clust.items():
					# check if sequence is assigned to a cluster, add values to dict
					data_dict[dist_thresh][prot_id][int(clust)][seq] = seq_2_values[seq]

	return data_dict, thresh_dict

def create_network( dist_idx, thresh_dict, ids, data_dict, outDf, vis_seed, net_out_dir, make_net_vis):
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
		cluster_sizes.append(len(data_dict[thresh_dict[id_][dist_idx]][id_][num]) * 25)

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
		weights = nx.get_edge_attributes(G,'weight').values()
		pos=nx.spring_layout(G, k=0.75, seed=vis_seed)
		fig, ax = plt.subplots(figsize = (20, 14))
		nx.draw_networkx_nodes(G, pos, ax=ax, node_size=cluster_sizes, node_color=cluster_colors)
		nx.draw_networkx_edges(G, pos, ax=ax, connectionstyle=f'arc3, rad = 0.25', arrows=True, width=list(weights))
		nx.draw_networkx_labels(G, pos, ax=ax)
		fig.savefig(f"{net_out_dir}/network_{dist_idx +1}_visualization.png", bbox_inches='tight', dpi=300)

	#---------Summary Statistics-----------
	# TODO: maybe put this in it's own function
	seq_data = list()

	orphaned_cluster_dict = defaultdict(set)
	for id_ in ids:
		for num in data_dict[thresh_dict[id_][dist_idx]][id_]:
			orphaned_cluster_dict[id_].add(num)

	counts_dict = defaultdict(list)
	net_stats = list()
	comp_dict = defaultdict(list)
	comps = nx.strongly_connected_components(G)
	for net_num, comp in enumerate(comps):
		# initialze dict for network num
		for id_ in ids:
			counts_dict[id_].append(0)
			comp_dict[id_].append(list())

		# find average linkage score
		avg_linkage_score = 0
		edge_count = 0

		# find number of sequences
		seq_num = 0

		for i, clust_1 in enumerate(comp):
			clust_1_list = clust_1.split("_")
			id_1 = clust_1_list[0]
			num_1 = int(clust_1_list[1])

			seq_num += len(data_dict[thresh_dict[id_1][dist_idx]][id_1][num_1])

			# assighn network number to sequence
			for seq in data_dict[thresh_dict[id_1][dist_idx]][id_1][num_1].keys():
				seq_data.append( (seq, net_num) )

			orphaned_cluster_dict[id_1].remove(num_1)

			counts_dict[id_1][net_num] += 1

			comp_dict[id_1][net_num].append(num_1)

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

		net_stats.append( (net_num, avg_linkage_score, seq_num) )

		for id_ in comp_dict.keys():
			comp_dict[id_][net_num] = ",".join(map(str, comp_dict[id_][net_num]))
	
	# get orphaned clusters stats
	orphaned_data = defaultdict()
	# go through each orphaned cluster
	for orp_id in orphaned_cluster_dict.keys():
		for orp_num in orphaned_cluster_dict[orp_id]:
			seq_num = 0
			net_num += 1
			# add a row
			for id_ in ids:
				counts_dict[id_].append(0)
				comp_dict[id_].append("")
			# increment orphaned cluster
			counts_dict[orp_id][net_num] += 1
			comp_dict[orp_id][net_num] = str(orp_num)

			for seq in data_dict[thresh_dict[orp_id][dist_idx]][orp_id][orp_num].keys():
				seq_data.append( (seq, net_num) )

			# app row to net stats
			net_stats.append( (net_num, 0, len(data_dict[thresh_dict[orp_id][dist_idx]][orp_id][orp_num])) )


	countsDf = pd.DataFrame.from_dict(counts_dict).rename_axis("MultiProteinCluster")
	compsDf = pd.DataFrame.from_dict(comp_dict).rename_axis("MultiProteinCluster")
	statsDf = pd.DataFrame( net_stats, columns=["MultiProteinCluster", "AvgNormLinkageScore", "NumberOfSequences"]).set_index("MultiProteinCluster")
	statsDf = countsDf.merge(statsDf, on="MultiProteinCluster", how="right")
	for id_ in ids:
		statsDf.rename( columns = {id_:f"{id_}_Cluster_Count"}, inplace = True)
	statsDf = statsDf.merge(compsDf, on="MultiProteinCluster", how="left")
	statsDf.to_csv(f"{net_out_dir}/network_{dist_idx +1}_summary_stats.tsv", sep="\t")

	seqDf = pd.DataFrame(seq_data, columns=["Sequence", "MultiProteinCluster"]).set_index("Sequence")
	seqDf.to_csv(f"{net_out_dir}/network_{dist_idx +1}_multiprotein_cluster_sequences.tsv", sep="\t")
	#--------------------------------------

#------------------------------------------
if __name__ == "__main__":
	main();