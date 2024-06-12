#!/usr/bin/env python

from generateClusters import cluster
from clusterLinkage import find_linkage_scores
import argparse
import os
import glob
import tempfile
import pandas as pd
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(description="Generates clusters and creates linkage network", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	# Required arguments
	parser.add_argument("-i", "--input-dir", type=str, required=True, help="Directory containing protein directories with fasta files in each.")
	parser.add_argument("-m", "--metadata", type=str, required=True, help="Tab-delimited metadata file that links sequences to info that can be used to link sequences across clusters.")
	parser.add_argument("-d", "--distance-thresh", nargs="+", type=float, required=False, help="Required for Hierarchical Clustering. Distance thresholds to use for hierarchical clustering. Multiple values may be provided, all of which should be between 0 and 1.")
	
	# Optional arguments (generate clusters)
	parser.add_argument("--method", type=str, required=False, default="Hierarchical", help="Method of generating clusters. Current options are: Hierarchical or Louvain")
	parser.add_argument("-k", "--kmer-size", type=int, required=False, default=7, help="Size of kmers used to compare sequences.")
	parser.add_argument("-p", "--min-propn", type=float, required=False, default=0, help="Proportion of the top 10%% of sequence sizes to be included in the initial round of clustering.")
	parser.add_argument("--make-dist-boxplots", required=False, default=False, action="store_true", help="Optional tab-delimited file that can be used to link the input sequences to metadata. If provided, summary statistics about the generated clusters will be generated.")
	parser.add_argument("--min-seqs", type=int, default=2, required=False, help="Number of sequences a cluster needs to be included in the boxplot.")
	parser.add_argument("--make-sim-hists", required=False, default=False, action="store_true", help="If provided, distribution of k-mer similarities for each input file will be provided.")
	parser.add_argument("--generate-vis", required=False, default=False, action="store_true", help="If provided, generates visualization of clusters for clustering method.")
	
	# Optional arguments (cluster linkage)
	parser.add_argument("-t", "--thresh-matrix", type=str, required=False, default="", help="Filepath to tab delimited file of the threshold matrix. Each column should represent the same protein id and the transposed rows "
																							"of the manifest file. Each row should be each seperate linkage network to create. The entries are the thresholds to use for each "
																							"protein id for each network to create. Optionally, an entry can be left blank and it will not be contributed to that network")
	parser.add_argument("--seq-name-header", default="SequenceName", type=str, required=False, help="Header for sequence name column in the metadata file.")
	parser.add_argument("--linkage-cols", nargs="+", default=["NCBIaccession-NT"], type=str, required=False, help="Header of columns in the metadata file t consider for linkage. "
																											"Any given sequence only contributes one point to a linkage score")
	parser.add_argument("--col-val-delim", default=",", type=str, required=False, help="Delimiter for multiple values of a column associated with a single sequence in the metadata file.")
	parser.add_argument("--make-network-vis", required=False, default=False, action="store_true", help="If provided, network visualization will be created for each threshold. The size of the nodes are based on the "
																										"number of sequences in the cluster and the size of the edges are based on the normalized linkage score.")
	parser.add_argument("--vis-seed", type=int, required=False, default="1234", help="Seed used to generate network visualization. Changes node positions.")

	args=parser.parse_args()

	if not os.path.exists(args.input_dir):
		raise Exception("Input directory does not exist")

	net_species_dict = defaultdict(lambda: pd.DataFrame())

	method_cluster_name = f"{args.method.lower()}_clusters"
	multi_prot_net_name = "multiprotein_networks"

	for folder in os.listdir(args.input_dir):
		prot_dir = os.path.join(args.input_dir, folder)
		if os.path.isdir(prot_dir) and len(glob.glob(os.path.join(prot_dir, '*.fasta'))) > 0:
			gen_clust_input = list()

			for file in glob.glob(os.path.join(prot_dir, '*.fasta')):
				gen_clust_input.append(file)

			cluster(
				input_files = gen_clust_input,
				method = args.method,
				distance_thresh = args.distance_thresh,
				kmer_size = args.kmer_size,
				min_propn = args.min_propn,
				make_dist_boxplots = args.make_dist_boxplots,
				boxplots_output_dir = os.path.join(prot_dir, f"{method_cluster_name}_boxplots"),
				min_seqs = args.min_seqs,
				make_sim_hists = args.make_sim_hists,
				hists_output_dir = os.path.join(prot_dir, f"{method_cluster_name}_histograms"),
				gen_vis = args.generate_vis,
				vis_output_dir = os.path.join(prot_dir, f"{method_cluster_name}_visualization"),
				output_dir = os.path.join(prot_dir, method_cluster_name)
				)

			# create temp manifest file
			id_dict = dict()
			for file in glob.glob(os.path.join(prot_dir, method_cluster_name, "*.fasta.tsv")):
				id_dict[file[len(os.path.join(prot_dir, method_cluster_name, f"clusters_{folder}_")) : len(file) - len(".fasta.tsv")]] = file

			with tempfile.NamedTemporaryFile( suffix=".tsv" ) as tempManifest:
				pd.DataFrame(id_dict.items(), columns=["ProtID", "ClustersFile"]).to_csv(tempManifest.name, sep="\t", index=False)

				find_linkage_scores(
					manifest_file = tempManifest.name,
					metadata = args.metadata,
					thresh_matrix = args.thresh_matrix,
					seq_header = args.seq_name_header,
					linkage_cols = args.linkage_cols,
					col_val_delim = args.col_val_delim,
					make_net_vis = args.make_network_vis,
					vis_seed = args.vis_seed,
					output_dir = os.path.join(prot_dir, multi_prot_net_name)
					)

			# append to network species
			for net in glob.glob(os.path.join(prot_dir, multi_prot_net_name, "*")):
				if os.path.isdir(net):
					network_name = os.path.basename(net)
					curr_spec_df = pd.read_csv(os.path.join(net, f"{network_name}_multiprotein_cluster_sequences.tsv"), sep="\t")
					curr_spec_df["MultiProteinCluster"] = f"{folder}_" + curr_spec_df["MultiProteinCluster"].astype(str)
					net_species_dict[network_name] = pd.concat([net_species_dict[network_name], curr_spec_df], ignore_index=True)

	# save networks
	merged_spec_dir = os.path.join(args.input_dir, "merged_sequences")
	if not os.path.exists(merged_spec_dir):
		os.mkdir(merged_spec_dir)

	for net in net_species_dict.keys():
		net_species_dict[net].to_csv(os.path.join(merged_spec_dir, f"{net}_merged_multiprotein_cluster_sequences.tsv"), sep="\t", index=False)


if __name__ == "__main__":
	main()

