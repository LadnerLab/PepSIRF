#!/usr/bin/env python

import optparse
import itertools as it
import numpy as np
import matplotlib.pyplot as plt


# This script can be used to generate count matrices that normalize using negative controls

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--input',  help='Matrix of normalized (e.g., column-sum) read counts. [None]')
    p.add_option('--negMat',  help='Optionally, this flag can be used to supply a separate data matrix containing SB samples. [None]')
    p.add_option('-s', '--start',  help='Optional approach for identifying negative controls. Can provide unique string at the start of all negative control samples [None]')
    p.add_option('-n', '--negNames', help='Optional approach for identifying negative controls. Comma-separated list of negative control samples names [None]')
    p.add_option('-a', '--approach', help='Approach for normalization. Options: "ratio" and "diff" [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file, which will contain ratios of sample Norm Count/average negative control sample Norm Count [None, REQ]')
    opts, args = p.parse_args()
    
    inD = parseCounts(opts.input)
    
    if opts.negMat:
        sbD = parseCounts(opts.negMat)
    
    if opts.start and opts.negNames:
        print("Script does not currently support a combination of '--start' and '--negNames', using '--negNames'")
        negL = opts.negNames.split(",")
    elif opts.negNames:
        negL = opts.negNames.split(",")
    elif opts.start:
        if opts.negMat:
            negL = [x for x in sbD if x.startswith(opts.start)]
        else:
            negL = [x for x in inD if x.startswith(opts.start)]
    else:
        negL = False
        print("Must use either '--start' or '--negNames' to specify the negative control samples to use for normalization")
    
    if negL:
        samps = list(inD.keys())
        negNormD = {x:{} for x in samps}
        for p in inD[samps[0]]:
            if opts.negMat:
                avgNeg = np.mean([sbD[x][p] for x in negL])
            else:
                avgNeg = np.mean([inD[x][p] for x in negL])
            
            for s in samps:
                if opts.approach == "ratio":
                    negNormD[s][p] = inD[s][p]/avgNeg
                elif opts.approach == "diff":
                    negNormD[s][p] = inD[s][p] - avgNeg
                else:
                    print("%s is not a recognized approach" % opts.approach)
                    
    writeCounts(negNormD, opts.out)
    
#----------------------End of main()

def writeCounts(cd, outname):
    probeNames = sorted(cd[list(cd.keys())[0]].keys())
    sampNames =  sorted(list(cd.keys()))
    with open(outname, "w") as fout:
        fout.write("Probe\t%s\n" % ("\t".join(sampNames)))
        for each in probeNames:
            fout.write("%s\t%s\n" % (each, "\t".join([str(cd[x][each]) for x in sampNames])))

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

