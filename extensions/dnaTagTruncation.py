#!/usr/bin/env python

import argparse, os
import inout as io
import fastatools as ft
from collections import defaultdict    

# The purpose of this script is to create truncated versions of the DNA tags used as a reference for demux
# The output will include the trncated reference file and a file containing the tag names that are still unique after the truncation

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")

    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-l', '--lens', help='Desired length(s) for truncated sequences. Can be a single integer, or multiple integers separated by commas', required=True)
    reqArgs.add_argument('-r', '--refs',  help='Input reference sequences to be truncated, in fasta format.', required=True)

    args = parser.parse_args()

    #Read in target fasta file
    names, seqs = ft.read_fasta_lists(args.refs)
    
    #Save base file name, from which output file names will be derived
    baseName = os.path.basename(args.refs)
    
    #Parse desired truncation length(s)
    targetLens = [int(x) for x in args.lens.split(",")]

    # Step through each length and create output files
    for each in targetLens:
        
        #Generate truncated sequences
        newSeqs = [s[:each] for s in seqs]
        
        #Write out new fasta file
        ft.write_fasta(names, newSeqs, "trunc%d_%s" % (each, baseName))
        
        #Count the number of times each sequence occurs in the list
        countD = defaultdict(list)
        for i, s in enumerate(newSeqs):
            countD[s].append(names[i])
        
        # Write out file with non-unique names, one row per sequence that is non-unique
        with open("trunc%d_%s_notUniqGrouped.tsv" % (each, baseName.split(".")[0]), "w") as foutG:
            for k,v in countD.items():
                if len(v) >1:
                    foutG.write("%s\t%s\n" % (k, "\t".join(v)))
        
        #Make list of sequences that are still unique, and another for those that aren't
        uniq = []
        notUniq = []
        for s,nl in countD.items():
            if len(nl)==1:
                 uniq.append(nl[0])
            elif len(nl)>1:
                notUniq+=nl
            else:
                print("Something went wrong. The name list for %s is empty." % (s))
        
        #Write out list of names for unique seqs
        io.writeList(uniq, "trunc%d_%s_uniq.txt" % (each, baseName.split(".")[0]))

        #Write out list of names for non-unique seqs
        io.writeList(notUniq, "trunc%d_%s_notUniq.txt" % (each, baseName.split(".")[0]))
        

#----------------------End of main()

    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

