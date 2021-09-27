#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

from collections import defaultdict

import optparse, sys, random
import matrixtools as mt             #Available at https://github.com/jtladner/Modules
#import itertools as it
#import inout as io             #Available at https://github.com/jtladner/Modules
#import numpy as np

# This script is used for downsampling a PepSeq dataset

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-d', '--data',  help='Matrix containing reads counts to downsample. [None, REQ]')
    p.add_option('-s', '--sampleName',  help='Name of focal sample contained within matrix. This is the sample that will serve as the starting place [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")
    p.add_option('-o', '--out',  help='Name for output file containing subset datasets. Output names are standardized and NOT related to the starting sample name. [None, REQ]')
    p.add_option('--min', type="int", help="Minimum downsampled dataset size. [None, REQ]")
    p.add_option('--max', type="int", help="Maximum downsampled dataset size. [None, REQ]")
    p.add_option('--step', type="int", help="Step size between dataset sizes. [None, REQ]")
    p.add_option('--reps', type="int", help="Number of replicates to generate at each size. [None, REQ]")

    opts, args = p.parse_args()
        
    # Print command used 
    print("Command run: '%s'" % ("  ".join(sys.argv)))
    
    #Read in data file 
    dataD = mt.parseCounts(opts.data, delim=opts.delim)
    
    # Check to make sure that the provided sample name is present in the provided data matrix
    if opts.sampleName in dataD:
        #Generate probability dictionary for focal sample
        probDict = countProbDict(dataD[opts.sampleName])
        
        #Initiate dictionary to hold downsampled data, also include complete dataset as "Full"
        toOut = {"Full":dataD[opts.sampleName]}
        
        #Generate downsampled datasets
        depths = list(range(opts.min, opts.max + 1, opts.step))
        for d in depths:
            for i in range(opts.reps):
                thisD = {k:0 for k in toOut["Full"]}
                sampleCounts(probDict, d, thisD)
                toOut["%d_%d" % (d, i)] = thisD
        
        mt.writeCounts(toOut, opts.out)
        
    else:
        print("%s NOT found in provided dataset!" % (opts.sampleName))
#----------------------End of main()

def countProbDict(infoD):
    num = 0
    cpD = {}
    for each, count in infoD.items():
        for i in range(int(count)):
            cpD[num] = each
            num+=1
    return cpD

def sampleCounts(probD, num2samp, d2fill):
    keys2samp = random.sample(range(len(probD)), num2samp)
    for k in keys2samp:
        d2fill[probD[k]]+=1


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

