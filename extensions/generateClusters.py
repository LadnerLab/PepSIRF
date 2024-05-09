import scipy.cluster.hierarchy as sch
import pandas as pd
import numpy as np
import kmertools as kt     # Available from https://github.com/jtladner/Modules
import inout as io         # Available from https://github.com/jtladner/Modules

import argparse
import os
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(description="Generates clusters of similar sequences for each protein")
	parser.add_argument("--meta-filepath", type=str, metavar="", required=True, help="")
	parser.add_argument("--input-files", nargs="+", default=[], metavar="", required=True, help="")
	parser.add_argument("--distance-thresh", nargs="+", type=float, metavar="", required=True, help="")
	parser.add_argument("--kmer-size", type=int, metavar="", required=True, help="")
	parser.add_argument("--min-propn", type=float, metavar="", required=False, default=1, help="Proportion of the max sequence size to not cluster.")
	parser.add_argument("-o", "--output-dir", type=str, metavar="", required=True, help="Directory to save cluster files.")
	args=parser.parse_args()

	cluster(
		meta_filepath = args.meta_filepath,
		input_files = args.input_files,
		distance_thresh = args.distance_thresh,
		kmer_size = args.kmer_size,
		min_propn = args.min_propn,
		output_dir = args.output_dir
		)

def cluster(
	meta_filepath: str,
	input_files: list,
	distance_thresh: list,
	kmer_size: int,
	min_propn: float,
	output_dir
	) -> None:
	
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	else:
		print(f"Warning: The directory \"{output_dir}\" already exists. Files may be overwritten.")

	# meta_filepath = "PM1_targets_taxInfo.tsv"
	speciesD = io.fileDictHeader(meta_filepath, "SequenceName", "Species")
	sidD = io.fileDictHeader(meta_filepath, "SequenceName", "SpeciesID")
	sid2spD = io.fileDictHeader(meta_filepath, "SpeciesID", "Species")

	'''
	# I want the user to be able to provide multiple input files, each which will be clustered separately
	input_files = ["39733_Astroviridae_Targets.fasta", "2842321_Kolmioviridae_Targets.fasta"]
	# I want the user to be able to provide multiple distance thresholds to use for clustering
	distance_thresh = [0.8, 0.9, 0.95]
	# I want the user to specify the kmer size to use
	kmer_size=7
	# I want the user to specify the location of an output directory
	# If this directory doesn't already exist, please create it
	output_dir = "clusters"
	'''

	for i in input_files:
	    print(i) # This doesn't need to be done in the standalone script, just temp for tracking
	    outD = {}


	    kmer_dict = kt.kmerDictSetFasta(i,kmer_size,["X"])
	    dists, seqNames = calcDistances(kmer_dict)
	    
	    for dt in distance_thresh:
	        clusters = clusterSeqs(dists,dt, seqNames)
	        
	        #Make a dictionary linking a sequence name to its assigned cluster
	        seq2clust = {}
	        for clustNum, seqL in clusters.items():
	            for sn in seqL:
	                seq2clust[sn] = clustNum
	        #Add cluster numbers to putput dictionary in order of seqNames
	        outD[dt] = [seq2clust[sn] for sn in seqNames]
	                
	        # This is an example of calculating some summary statistics to help us better understand the clusters
	        # We will likely want to expand this functionality in the future
	        numClusts = len(clusters)
	        multi, initialSpecies, multiClustSpecies = clustersBySpecies(clusters, speciesD)
	        print(f"{dt}\t{numClusts}\t{multi}\t{multi/numClusts:.4f}\t{len(multiClustSpecies)}\t{len(multiClustSpecies)/initialSpecies:.4f}")
	    
	    # Write out results for this input file
	    outDF = pd.DataFrame(outD, index=seqNames)
	    outDF.to_csv(f"{output_dir}/clusters_{os.path.basename(i)}.tsv", sep="\t", index_label="Sequence")

def calcDistances(kD):
    seqNames = list(kD.keys())

    # Currently, one distance measure is implemented
    # But, in the future, it would be nice to support multiple distance options
    dists = []
    for i, ni in enumerate(seqNames):
        ki = kD[seqNames[i]]
        for j, nj in enumerate(seqNames):
            if j>i:
                kj = kD[seqNames[j]]
                ovlp = ki.intersection(kj)
                dists.append(1-(len(ovlp)/(min([len(ki), len(kj)]))))
    return dists, seqNames

def clusterSeqs(dists, distThresh, seqNames, meth='average'):
    
    hm = sch.linkage(np.array(dists), method=meth)
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
#            print(spS)
        for each in spS:
            clusterCountD[each]+=1
    
    return multi, len(clusterCountD),{k:v for k,v in clusterCountD.items() if v>1}

if __name__ == "__main__":
	main()