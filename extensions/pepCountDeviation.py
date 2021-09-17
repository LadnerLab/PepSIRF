#!/usr/bin/env python

#import matplotlib.pyplot as plt
from scipy.stats import poisson
import optparse, os
#import numpy as np
#import inout as io
from collections import defaultdict

# This script reads in normalized read counts and zscores and generates a plot comparing them

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-r', '--raw',  help='A matrix containing raw scores for samples of interest. [None, REQ]')
    p.add_option('-n', '--normControl',  help='A matrix containing normalized scores for the negative control samples, which will be used to calculate peptide relative abundance. [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file, which will contain an estimate of number of missing peptides, given read count, for each sample [None, REQ]')
    p.add_option('-u', '--upperValue', type='int', default=0, help='Upper value to consider in calculations. Values equal to or less than this value will be compared with expectations [0]')
    opts, args = p.parse_args()
    
    # Generate relative probabilities of different peptides
    normD = parseCounts(opts.normControl)
    negSamps = list(normD.keys())
    probD = probeProb(normD, negSamps)
    numPeps = len(normD[negSamps[0]])
    
    print("SB Probs Done!")
    
    # Read in raw counts
    rawD = parseCounts(opts.raw)

    print("Read in raw counts!")
    
    counter=0
    with open("%s" % opts.out, "w") as fout:
        fout.write("Sample\tMissingPeptides\n")
        for s, inf in rawD.items():
            rc = sum(inf.values())
            zeros = sum([1 for x in inf.values() if x<=opts.upperValue])
            missing = expectUniq(probD, rc, opts) - (numPeps-zeros)
            fout.write("%s\t%d\n" % (s,missing))
            counter+=1
            print(counter, s, zeros, missing)

#----------------------End of main()


def expectUniq(probD, readCount, opts):
    totProb =  sum([1-poisson.cdf(opts.upperValue, readCount*(v)) for k,v in probD.items()])
    return totProb

def probeProb(normD, controls):
    sumD={}
    total=0
    for pep in normD[controls[0]]:
        thisSum = sum([normD[x][pep] for x in controls])
        total+=thisSum
        sumD[pep] = thisSum
    return {k:(v/total) for k,v in sumD.items()}

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

