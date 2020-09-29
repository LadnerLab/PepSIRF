#!/usr/bin/env python

import optparse, os
import numpy as np
from collections import defaultdict

# This script reads in normalized read counts and zscores and generates a plot comparing them

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-s', '--samples',  help='The name of the file containing sample pairs, denoting which samples are in duplicate. This file must be tab-delimited with one pair or samples per line. [None, REQ]')
    p.add_option('-t', '--thresh',  help='Tab-delimited file with information on score matrices and thresholds to use for determining enriched peptides. One row per matrix, 1st column = score matrix file, 2nd column = threshold(s), comma-sep if multiple [None, REQ]')
    p.add_option('-r', '--raw_scores',  help='Optionally, a raw count matrix can be included. This matrix must contain the raw counts of each probe. If included, "min_raw_count" must also be specified. [None, REQ]')
    p.add_option('--rc', '--raw_score_constraint', type="float", help="The minimum raw count a sample can have for all of its peptides in order for any of the probes in that sample to be considered enriched. The sum of each probe's raw count in each sample must be at least this value in order for the sample to be considered.  [None]")
    p.add_option('-d', '--outDir', help='Directory name for output files. Will be created. [None]')
    p.add_option('-o', '--outfile_suffix', help='Suffix for output files [None, REQ]')
    opts, args = p.parse_args()
    
    # Read in matrices and thresholds
    matDL = []
    threshLL = []
    with open(opts.thresh, "r") as fin:
        for line in fin:
            matrixF,thisThresh = line.rstrip("\n").split("\t")
            matDL.append(parseCounts(matrixF))
            threshLL.append([float(x) for x in thisThresh.split(",")])

    if opts.raw_scores:
        rD = parseCounts(opts.raw_scores)
        tooFew = [samp for samp,inf in rD.items() if sum(inf.values())<opts.rc]
        for each in tooFew:
            for mD in matDL:
                del(mD[each])
    
    # Read in pairs
    pairs = []
    with open(opts.samples, "r") as fin:
        for line in fin:
            pairs.append(line.rstrip("\n").split("\t"))
    
    #Create output dir
    if os.path.isdir(opts.outDir):
        print("The output directory %s already exists!" % (opts.outDir))
    else:
        os.mkdir(opts.outDir)
        
    #Generate enricehd probe lists
    outLists = defaultdict(list)
    for a,b in pairs:
        if a in matDL[0] and b in matDL[0]:
            for p in matDL[0][a]:
                vals = [[D[a][p], D[b][p]] for D in matDL]
                if sum([meetsThresh(vals[i], threshLL[i]) for i in range(len(vals))]) == len(vals):
                    outLists[(a,b)].append(p)
        else:
            print("At least one in this pair is missing from matrix, or has too few total reads: %s, %s" % (a,b))

    #Write output files
    for pair, probes in outLists.items():
        if len(probes)>1:
            with open("%s/%s%s" % (opts.outDir, "~".join(pair), opts.outfile_suffix), "w") as fout:
                fout.write("%s\n" % ("\n".join(probes)))

#----------------------End of main()

# For a single criteria, check that thresholds have been met
def meetsThresh(dataL, threshL):
    if max(dataL) >= max(threshL) and min(dataL) >= min(threshL):
        return True
    else:
        return False


def parseCounts(countFile, delim="\t"):
    counts={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc == 1:
                names = cols[1:]
                for n in names:
                    counts[n]={}
            else:
                for i, count in enumerate(cols[1:]):
                    counts[names[i]][cols[0]] = float(count)
    return counts


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

