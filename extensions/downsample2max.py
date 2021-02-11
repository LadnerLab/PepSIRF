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

# This script is used for to downsample a a full raw count PepSeq matrix
# The user provides a max number of reads to include for each sample
# The user can also specify a minimum read count for inclusion in the output matrix

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-d', '--data',  help='Matrix containing reads counts to downsample. [None, REQ]')
    p.add_option('-o', '--out',  help='Name for output file containing subset datasets. Output names are standardized and NOT related to the starting sample name. [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")
    p.add_option('--max', type="int", help="Size of output downsampled datasets. [None, REQ]")
    p.add_option('--min', type="int", default=0, help="Minimum dataset size for inclusion in the output. Samples with dataset sizes between min and max will be included without downsampling [None, REQ]")

    opts, args = p.parse_args()
        
    # Print command used 
    print("Command run: '%s'" % ("  ".join(sys.argv)))
    
    #Read in data file 
    dataD = mt.parseCounts(opts.data, delim=opts.delim)
    
    #Initiate dictionary to hold downsampled data
    toOut = {}
    
    # Step through each sample in the input data matrix
    for s, dd in dataD.items():
        
        numReads = sum(dd.values())
        
        #If the total number of reads is greater than the max set by user, then downsample
        if numReads > opts.max:

            #Generate probability dictionary for focal sample
            probDict = countProbDict(dd)
    
            #Generate downsampled dataset
            thisD = {k:0 for k in dd}
            sampleCounts(probDict, opts.max, thisD)
            toOut[s] = thisD
    
        elif numReads >= opts.min:
            toOut[s] = dd
    
    mt.writeCounts(toOut, opts.out)
    
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

