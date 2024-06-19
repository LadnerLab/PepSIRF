#!/usr/bin/env python

import argparse, os
import fastatools as ft # Available from https://github.com/jtladner/Modules
import pandas as pd
import inout as io         # Available from https://github.com/jtladner/Modules
from collections import defaultdict

#This script is designed to generate subset fasta files to serve as inputs for sequence-based clustering


def main():
    parser = argparse.ArgumentParser(description="Generates subset fasta files based on metadata")
    
    # Required arguments
    parser.add_argument("-i", "--input", metavar="", required=True, help="Fasta file containing sequences of interest.")
    parser.add_argument("-s", "--subsets", metavar="", required=True, help="Tab-delimited file with metadata for identifying subsets of interest.")
    parser.add_argument("-m", "--meta-filepath", type=str, metavar="", required=True, help="Tab-delimited file that can be used to link the input sequences to metadata.")

    args=parser.parse_args()

    #Read in fasta file as dictionary
    fD = ft.read_fasta_dict(args.input)
    
    #Read in sequence metadata file as a dataframe
    mDF = pd.read_csv(args.meta_filepath, sep="\t", header=0, index_col=0, keep_default_na=False)
    
    # Indexes for sequences with missing families
    missFamIdx = mDF[mDF["Family"]==""].index
    # Replace all of these with Unclassified
    for idx in missFamIdx:
        mDF.at[idx, 'Family'] = "Unclassified"

    #Read in cluster metadata
    cD = {"Family":defaultdict(dict), "Genus":defaultdict(dict)}
    with open (args.subsets, "r") as fin:
        next(fin)
        for line in fin:
            cols=line.rstrip("\n").split("\t")
            cD[cols[0]][cols[1]][cols[2]] = cols[3]

    # Step through each family of interest
    for each, subCD in cD["Family"].items():
        io.makeDir(each)
        famDF = mDF[mDF["Family"]==each]
        
        #Dictionary to link sequence names to cluster groups
        seqD = defaultdict(list)
        
        #Link sequence names to cluster groups
        for i, row in famDF.iterrows():
            for eachP in row["Protein"].split(","):
                seqD[cD["Family"][each][eachP]].append(i)
    
        #Write out a fasta file for each cluster group
        for clustGroup, seqL in seqD.items():
            # To handle duplicates
            seqL = sorted(list(set(seqL)))
            ft.write_fasta(seqL, [fD[sn] for sn in seqL], f"{each}/{each}_{clustGroup}.fasta")

    # Step through each genus of interest
    for each, subCD in cD["Genus"].items():
        io.makeDir(each)
        famDF = mDF[mDF["Genus"]==each]
        
        #Dictionary to link sequence names to cluster groups
        seqD = defaultdict(list)
        
        #Link sequence names to cluster groups
        for i, row in famDF.iterrows():
            for eachP in row["Protein"].split(","):
                seqD[cD["Genus"][each][eachP]].append(i)
    
        #Write out a fasta file for each cluster group
        for clustGroup, seqL in seqD.items():
            # To handle duplicates
            seqL = sorted(list(set(seqL)))
            ft.write_fasta(seqL, [fD[sn] for sn in seqL], f"{each}/{each}_{clustGroup}.fasta")



###---------------End of main()-----------------------------------


if __name__ == "__main__":
    main()