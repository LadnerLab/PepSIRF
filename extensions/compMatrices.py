#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import optparse, glob, os, sys
#import inout as io             #Available at https://github.com/jtladner/Modules
import matrixtools as mt             #Available at https://github.com/jtladner/Modules
import numpy as np

from collections import defaultdict

# This script is used to compare two matrices to make sure they are identical or near-identical

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    #Experimental data
    p.add_option('-1', '--m1',  help='First matrix [None, REQ]')
    p.add_option('-2', '--m2',  help='Second matrix [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")
#    p.add_option('--minPosSlope', type="float", default=100, help="If the slope for the positive control sample is less than this, a warning will be issued. [100]")
#    p.add_option('--maxNegSlope', type="float", default=30, help="If the slope for the negative control sample is greater than this, a warning will be issued. [30]")
#    p.add_option('--minRvalue', type="float", default=0.97, help="If the correlation coefficient for a sample is less than this, a warning will be issued. [0.97]")
#    p.add_option('--maxSlopeRatio', type="float", default=1.2, help="If the ratio of the (slope from empty control and most dilute sample)/(best dilution slope) is greater than this, a warning will be issued. [1.2]")

    opts, args = p.parse_args()
        
    #Read in matrices
    m1D = mt.parseCounts(opts.m1, delim=opts.delim, valType="str")
    m2D = mt.parseCounts(opts.m2, delim=opts.delim, valType="str")
    
    keepTesting = checkAsStr(m1D, m2D)
    
    if keepTesting:
        diffs, diffprops = calcDiffs(m1D, m2D)
        print("Differences:\n\tMax:%.9f\n\tMin:%.9f\n\tAvg:%.9f" % (max(diffs), min(diffs), np.mean(diffs)))
        print("Difference Proportions:\n\tMax:%.9f\n\tMin:%.9f\n\tAvg:%.9f" % (max(diffprops), min(diffprops), np.mean(diffprops)))
#        testAsFloats = checkAsStrTruncate(m1D, m2D)

#----------------------End of main()

def calcDiffs(m1D, m2D):
    diffs=[]
    diffprops=[]
    
    nan1=0
    nan2=0
    nanBoth=0
    inf1=0
    inf2=0
    infBoth=0
    
    for samp in m1D:
        for p in m1D[samp]:
            strVals = [m1D[samp][p], m2D[samp][p]]
            if 'inf' in strVals or '-inf' in strVals or 'nan' in strVals or '-nan' in strVals:
                if strVals[0] == 'inf' or strVals[0] == '-inf':
                    if strVals[1] == 'inf' or strVals[1] == '-inf':
                        infBoth+=1
                    else:
                        inf1+=1
                elif strVals[1] == 'inf' or strVals[1] == '-inf':
                    inf2+=1
                    
                if strVals[0] == 'nan' or strVals[0] == '-nan':
                    if strVals[1] == 'nan' or strVals[1] == '-nan':
                        nanBoth+=1
                    else:
                        nan1+=1
                elif strVals[1] == 'nan' or strVals[1] == '-nan':
                    nan2+=1
            else:
                vals = [float(m1D[samp][p]), float(m2D[samp][p])]
                diffs.append(max(vals)-min(vals))
                if max(vals) != 0:
                    diffprops.append(diffs[-1]/abs(max(vals)))
    
    print("Both with nan: %d, Both with inf: %d" % (nanBoth, infBoth))
    print("Only 1st with nan: %d, Only 1st with inf: %d" % (nan1, inf1))
    print("Only 2nd with nan: %d, Only 2nd with inf: %d" % (nan2, inf2))
    
    return diffs, diffprops
    

def findMaxDec(m1D, m2D):
    currMax=0
    for d in [m1D, m2D]:
        for s, sD in d.items():
            for p, v in sD.items():
                decNum = len(v.split(".")[1])
                if decNum > currMax:
                    currMax = decNum
    return currMax
    

def checkAsStrTruncate(m1D, m2D):
    maxDecimals = findMaxDec(m1D, m2D)
    decimals=0
    while decimals<=maxDecimals:
        for samp in m1D:
            for p in m1D[samp]:
                trunc1 = truncate(m1D[samp][p], decimals)
                trunc2 = truncate(m2D[samp][p], decimals)
                
                if trunc1 != trunc2:
                    print(samp, p, m1D[samp][p], m2D[samp][p])
                    print("Matrices are NOT IDENTICAL, even when only considering %d decimal places" % (decimals))
                    return True
        print("Matrices ARE IDENTICAL when only considering %d decimal places" % (decimals))
        decimals+=1
    return False

def truncate(string, decimals):
    parts=string.split(".")
    return parts[0] + parts[1][:decimals]

def checkAsStr(m1D, m2D):
    k1 = sorted(list(m1D.keys()))
    k2 = sorted(list(m2D.keys()))
    
    if k1==k2:
        print("Samples contained ARE IDENTICAL!")
        
        for samp in m1D:
            for p in m1D[samp]:
                if m1D[samp][p] != m2D[samp][p]:
                    print("Matrices are NOT IDENTICAL!")
                    return True
        
        print("Matrices are IDENTICAL!")
        return False

    else:
        print("Samples contained are NOT the same. Ending comparison!")
        print("Samples only seen in m1: %s" % (set(k1).difference(set(k2))))
        print("Samples only seen in m2: %s" % (set(k2).difference(set(k1))))
        return False

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

