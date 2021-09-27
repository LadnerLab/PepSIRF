#!/usr/bin/env python

import argparse, glob, os, sys
#import inout as io             #Available at https://github.com/jtladner/Modules
import fastatools as ft        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script reads in an aligned fasta file
# And generates a fasta file containing a consensus sequence

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--consName", default="consensus", help="Name to use for the consensus sequence in the output file")
    parser.add_argument("--includeTerminalDashes", default=True, action="store_false", help="By default, terminal '-' characters will not be considered in consensus generation.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-o", "--out", help="Output file name", required=True)
    reqArgs.add_argument("-a", "--aligned", help="Aligned fasta file", required=True)

    args = parser.parse_args()

    # Read in aligned fasta file
    names, seqs = ft.read_fasta_lists(args.aligned)
    
    # If ignoring terminal dashes, convert those dashes to '~' characters
    seqs = replaceTerminalDashes(seqs)
    
    # Create consensus sequence
    cons = makeConsensus(seqs)
    
    # Write out fasta with consensus
    ft.write_fasta([args.consName], [cons], args.out)
    
#----------------------End of main()

def replaceTerminalDashes(seqs, repChar="~"):

    new_seqs=[]
    for s in seqs:
        beg_count=0
        end_count=0
        for b in s:
            if b=='-': beg_count+=1
            else: break
        for b in s[::-1]:
            if b=='-': end_count+=1
            else: break

        if beg_count==0 and end_count==0: new_seqs.append(s)
        elif beg_count==0: new_seqs.append(s[:-end_count] + repChar*end_count)
        elif end_count==0: new_seqs.append(repChar*beg_count + s[beg_count:])
        else: new_seqs.append(repChar*beg_count + s[beg_count:-end_count] + repChar*end_count)

    return new_seqs

def makeConsensus(Seqs):
    cons=[]
    for p in range(len(Seqs[0])):
        bases = [x[p] for x in Seqs]
        uniq = list(set(bases).difference(set(["~"])))
        if len(uniq) > 1:
            uniq = list(set(bases).difference(set(["~", "X"])))
        top=""
        topC=0
        for b in uniq:
            bC = bases.count(b)
            if bC>topC:
                top=b
                topC=bC
        #To avoid positions where most sequences have a deletion
        if top != "-":
            cons.append(top)
    return "".join(cons)

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

