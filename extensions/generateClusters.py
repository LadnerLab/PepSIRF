#!/usr/bin/env python

import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import kmertools as kt	 # Available from https://github.com/jtladner/Modules
import fastatools as ft # Available from https://github.com/jtladner/Modules
import inout as io		 # Available from https://github.com/jtladner/Modules

import argparse, random, os
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(description="Generates clusters of similar sequences for each protein")
	
	# Required arguments
	parser.add_argument("-i", "--input-files", nargs="+", default=[], metavar="", required=True, help="Fasta files containing sequences to be clustered. One or more fasta files can be provided. Sequences in each will be separately clustered.")
	parser.add_argument("-d", "--distance-thresh", nargs="+", type=float, metavar="", required=True, help="Distance thresholds to use for hierarchical clustering. Multiple values may be provided, all of which should be between 0 and 1.")
	parser.add_argument("-o", "--output-dir", type=str, metavar="", required=True, help="Directory to save cluster files. This directory will be created, if it doesn't already exist.")

	# Optional arguments
	# parser.add_argument("--method", type=str, metavar="", required=False, default="Hierarchical", help="Method of generating clusters. Current options are: Hierarchical or Louvain")
	parser.add_argument("-k", "--kmer-size", type=int, metavar="", required=False, default=7, help="Size of kmers used to compare sequences.")
	parser.add_argument("-p", "--min-propn", type=float, metavar="", required=False, default=0, help="Proportion of the top 10%% of sequence sizes to be included in the initial round of clustering.")
	parser.add_argument("-m", "--meta-filepath", type=str, metavar="", required=False, help="Optional tab-delimited file that can be used to link the input sequences to metadata. If provided, summary statistics about the generated clusters will be generated.")
	parser.add_argument("--make-dist-boxplots", required=False, default=False, action="store_true", help="Optional tab-delimited file that can be used to link the input sequences to metadata. If provided, summary statistics about the generated clusters will be generated.")
	parser.add_argument("--boxplots-output-dir", type=str, metavar="", default="boxplots", required=False, help="Directory to save boxplots if boxplots are made.")
	parser.add_argument("--min-seqs", type=int, metavar="", default=2, required=False, help="Number of sequences a cluster needs to be included in the boxplot.")
	parser.add_argument("--make-sim-hists", required=False, default=False, action="store_true", help="If provided, distribution of k-mer similarities for each input file will be provided.")
	parser.add_argument("--hists-output-dir", type=str, metavar="", required=False, default="histograms", help="Directory to save histograms if histograms are made.")

	args=parser.parse_args()

	cluster(
		meta_filepath = args.meta_filepath,
		input_files = args.input_files,
		# method = args.method,
		distance_thresh = args.distance_thresh,
		kmer_size = args.kmer_size,
		min_propn = args.min_propn,
		make_dist_boxplots = args.make_dist_boxplots,
		boxplots_output_dir = args.boxplots_output_dir,
		min_seqs = args.min_seqs,
		make_sim_hists = args.make_sim_hists,
		hists_output_dir = args.hists_output_dir,
		output_dir = args.output_dir
		)

###---------------End of main()-----------------------------------


def cluster(
	meta_filepath: str,
	input_files: list,
	# method: str,
	distance_thresh: list,
	kmer_size: int,
	min_propn: float,
	make_dist_boxplots: bool,
	boxplots_output_dir: str,
	min_seqs: int,
	make_sim_hists: bool,
	hists_output_dir: str,
	output_dir: str
	) -> None:
	
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")
	if make_dist_boxplots:
		if not os.path.exists(boxplots_output_dir):
			os.mkdir(boxplots_output_dir)
		else:
			print(f"Warning: The directory \"{boxplots_output_dir}\" already exists. Boxplots may be overwritten.")
	if make_sim_hists:
		if not os.path.exists(hists_output_dir):
			os.mkdir(hists_output_dir)
		else:
			print(f"Warning: The directory \"{hists_output_dir}\" already exists. Histograms may be overwritten.")

# Want to make this section optional 
	# meta_filepath = "PM1_targets_taxInfo.tsv"
#	speciesD = io.fileDictHeader(meta_filepath, "SequenceName", "Species")
#	sidD = io.fileDictHeader(meta_filepath, "SequenceName", "SpeciesID")
#	sid2spD = io.fileDictHeader(meta_filepath, "SpeciesID", "Species")

	for i in input_files:
#		print(i) # This doesn't need to be done in the standalone script, just temp for tracking
		outD = {}
		
		# Read in sequences from fasta file
		fD = ft.read_fasta_dict(i)
		# Sort into separate dictionaries according to size threshold
		largeD, smallD = sortSeqsBySize(fD, min_propn)
		
		# Generate kmers for all sequences that meet size threshold
		kmer_dict = kt.kmerDictSet(largeD,kmer_size,["X"])
		# Calculate distances between all large sequences
		dists, similarities, seqNames = calcDistances(kmer_dict)

		# Generate kmers for all sequences that DO NOT meet size threshold
		kmer_dict_small = kt.kmerDictSet(smallD,kmer_size,["X"])
		shortSeqNames = list(smallD.keys())

		for dt in distance_thresh:
			# Cluster long sequences
			clusters = clusterSeqs(dists, dt, seqNames)
			
			if len(kmer_dict_small) > 0:
				# Start to assigning the short sequences to the clusters built using long sequences
				seq2clust = assignShortSequences(kmer_dict_small, kmer_dict, clusters, dt)
			else:
				# Make a dictionary linking a sequence name to its assigned cluster
				seq2clust = {}

			for clustNum, seqL in clusters.items():
				for sn in seqL:
					seq2clust[sn] = clustNum

			# Add cluster numbers to output dictionary in order of seqNames
			outD[dt] = [seq2clust[sn] for sn in seqNames + shortSeqNames]


			if make_dist_boxplots:
				make_boxplots(i, dt, clusters, kmer_dict, outD, boxplots_output_dir, min_seqs);
			if make_sim_hists:
				make_hist(i, similarities, hists_output_dir);


		# Write out results for this input file
		outDF = pd.DataFrame(outD, index=seqNames + shortSeqNames)
		outDF.to_csv(f"{output_dir}/clusters_{os.path.basename(i)}.tsv", sep="\t", index_label="Sequence")

		# Summary statistics
		clusterStats(outD, i , output_dir)

			# This is an example of calculating some summary statistics to help us better understand the clusters
			# We will likely want to expand this functionality in the future
# Want to make this section optional 
#			numClusts = len(clusters)
#			multi, initialSpecies, multiClustSpecies = clustersBySpecies(clusters, speciesD)
#			print(f"{dt}\t{numClusts}\t{multi}\t{multi/numClusts:.4f}\t{len(multiClustSpecies)}\t{len(multiClustSpecies)/initialSpecies:.4f}")

# For each short sequence, calculate the smallest distance between this sequence and the clusters of long sequences
	# Based on these distances, assign the short sequences to clusters
def assignShortSequences(shortKD, longKD, clusters, distThresh):
	clustAssign = {}
	for ssN, ssK in shortKD.items():
		distD = {}
		for clustNum, seqL in clusters.items():
			theseDists = []
			for lsN in seqL:
				lsK = longKD[lsN]
				ovlp = ssK.intersection(lsK)
				theseDists.append( 1-( len(ovlp)/(min([len(ssK), len(lsK)])) ) )
			distD[clustNum] = min(theseDists)
		
		bestScore = min(distD.values())
		if bestScore <= distThresh:
			bestClusts = [k for k,v in distD.items() if v==bestScore]
			clustAssign[ssN] = random.choice(bestClusts)
		else:
			clustAssign[ssN] = ""
	
	return clustAssign
		
def sortSeqsBySize(fastaDict, min_propn, topPerc=0.1):
	# Sort the sequence lengths in ascending order
	lenSortedL = sorted([len(seq) for seq in fastaDict.values()])
	# Determine the number of sequences needed to contain 10% of the total
	numForTopPerc = round(len(fastaDict)*topPerc)
	# Calculate the average size for the top specified percentage of sequences
	topAvg = np.mean(lenSortedL[-numForTopPerc:])
	# Calculate threshold for inclusion in the initial round of clustering
	thresh = topAvg*min_propn
	
	# Generate a dictionary to contain the longer sequences to be used in initial round of clustering
	largeD={}
	# Generate a dictionary to contain the shorter sequences that will not be used in initial round of clustering
	smallD={}
	for k,v in fastaDict.items():
		if len(v)>=thresh:
			largeD[k] = v
		else:
			smallD[k] = v
	
	return largeD, smallD

# within the same kmer dictionary
def calcDistances(kD):
	seqNames = list(kD.keys())

	# Currently, one distance measure is implemented
	# But, in the future, it would be nice to support multiple distance options
	dists = []
	similarities = []
	for i, ni in enumerate(seqNames):
		ki = kD[seqNames[i]]
		for j, nj in enumerate(seqNames):
			if j>i:
				kj = kD[seqNames[j]]
				ovlp = ki.intersection(kj)
				similarity = len(ovlp)/(min([len(ki), len(kj)]))
				similarities.append(similarity)
				dists.append(1-similarity)
	return dists, similarities, seqNames

#TODO: account for similarities
def calcDistancesBetweenClusters(kmer_dict1, kmer_dict2):
	seqNames1 = list(kmer_dict1.keys())
	seqNames2 = list(kmer_dict2.keys())

	dists = []
	for i, ni in enumerate(seqNames1):
		ki = kmer_dict1[seqNames1[i]]
		for j, nj in enumerate(seqNames2):
			kj = kmer_dict2[seqNames2[j]]
			ovlp = ki.intersection(kj)
			dists.append(1-(len(ovlp)/(min([len(ki), len(kj)]))))
	return dists


def clusterSeqs(dists, distThresh, seqNames, linkage_meth='average', clust_method):
	hm = sch.linkage(np.array(dists), method=linkage_meth)
	groups = sch.cut_tree(hm,height=distThresh)

	gD=defaultdict(list)
	for i,g in enumerate(groups):
		gD[g[0]].append(seqNames[i])

	return gD

def clustersBySpecies(clusters, speciesD):
	clusterCountD = defaultdict(int)
	multi=0

	for grpNum, memList in clusters.items():
		spL = [speciesD[m] for m in memList]
		spS = set(spL)
		if len(spS)>1:
			multi+=1
#			print(spS)
		for each in spS:
			clusterCountD[each]+=1
	
	return multi, len(clusterCountD),{k:v for k,v in clusterCountD.items() if v>1}

# Generate some summary statistics about the generated clusters
def clusterStats(outD, inName, output_dir):
	with open(f"{output_dir}/{inName}_summaryStats.tsv", "w") as fstats:
		fstats.write("Distance_Threshold\tNumberOfClusters\tAvgSeqsPerCluster\tSeqsUnassigned\n")
		for dt, clustL in outD.items():
			uniq = set(clustL).difference({""})
			numClust = len(uniq)
			numUnassigned = clustL.count("")
			perClust=np.mean([clustL.count(x) for x in uniq])
			fstats.write(f"{dt:.3f}\t{numClust}\t{perClust:.1f}\t{numUnassigned}\n")

def make_boxplots(input_filename, dist_thresh, clusters, kmer_dict, outD, boxplots_output_dir, min_seqs):
	clustDistWithin = list()
	clustDistBetween = list()

	for clustNum1, clustSeqNames1 in clusters.items():
		if len(clustSeqNames1) >= min_seqs:
			kmerSubDict1 = {seq:kmer_dict[seq] for seq in clustSeqNames1}

			clustDistsWithin, similarities, clustSeqNames = calcDistances(kmerSubDict1)

			for dist in clustDistsWithin:
				clustDistWithin.append((dist, clustNum1, "Within"))

			for clustNum2, clustSeqNames2 in clusters.items():
				if clustNum1 != clustNum2:
					kmerSubDict2 = {seq:kmer_dict[seq] for seq in clustSeqNames2}

					clustDistsBetween = calcDistancesBetweenClusters(kmerSubDict1, kmerSubDict2)

					for dist in clustDistsBetween:
						clustDistBetween.append((dist, clustNum1, "Between"))

	withinDf = pd.DataFrame(clustDistWithin, columns = ["Distance", "Cluster", "Type"])
	betweenDf=pd.DataFrame(clustDistBetween, columns = ["Distance", "Cluster", "Type"])

	clustDistDf = pd.concat([withinDf, betweenDf])

	width = max(outD[dist_thresh])
	height = 20
	fontsize = width * 5/4
	fig, ax = plt.subplots(figsize=(width, height), facecolor='w')
	sns.boxplot(x="Cluster", y="Distance", hue="Type", data=clustDistDf, ax=ax)
	ax.set_xlabel("Cluster", fontsize=fontsize)
	ax.set_ylabel("Distance", fontsize=fontsize)
	plt.grid()
	plt.savefig(f"{boxplots_output_dir}/{input_filename}_{dist_thresh}_boxplot.png", dpi=300, bbox_inches='tight')

def make_hist(input_filename, similarities, output_dir):
	fig, ax = plt.subplots(figsize=(11,8.5), facecolor='w')
	ax.hist(similarities, bins=100)
	plt.xticks(np.arange(0,1, step=0.1))
	plt.xlabel("similarity")
	plt.ylabel("count")
	plt.savefig(f"{output_dir}/{input_filename}_histogram.png")

if __name__ == "__main__":
	main()