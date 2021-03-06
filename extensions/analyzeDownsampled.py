#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

#from collections import defaultdict

import optparse, sys, random, glob
import matrixtools as mt             #Available at https://github.com/jtladner/Modules
import inout as io

# This script generates a set of plots to help with interpretation of downsampled datasets generated with "downsampleReads.py"

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-d', '--data',  help='Downsample data matrix generated by "downsampleReads.py". [None, REQ]')
    p.add_option('-e', '--enrDir',  help='Directory containing enriched peptide files, if available". [None, REQ]')
    p.add_option('-o', '--out',  help='Downsample data matrix generated by "downsampleReads.py". [None, REQ]')

    opts, args = p.parse_args()
        
    # Print command used 
    print("Command run: '%s'" % ("  ".join(sys.argv)))
    
    #Read in data file 
    dataD = mt.parseCounts(opts.data)
    
    #Cumulative norm diff figure
    fig,ax = plt.subplots(figsize=(6,5),facecolor='w')

    x=[]
    y=[]

    for s, d in dataD.items():
        if s != "Full":
            x.append(int(s.split("_")[0]))
            y.append(cumulativeNormDiff(dataD["Full"], d))

    ax.scatter(x, y, alpha=0.9, zorder=2, s=10)
    ax.set_xlabel("# Reads", fontsize=15)
    ax.set_ylabel("Cumulative Norm Difference", fontsize=15)
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    fig.savefig("%s_cumulativeNormDiff.png" % (opts.out), dpi=300, bbox_inches='tight')
    
    # Norm diff boxplots
    fig,ax = plt.subplots(figsize=(20,5),facecolor='w')

    raw=[]

    for s, d in dataD.items():
        if s != "Full" and s.split("_")[1] == "0":
            raw.append((int(s.split("_")[0]), [dataD["Full"][p] - d[p] for p in d]))

    y = [r[1] for r in sorted(raw)]
    
    ax.boxplot(y)
    ax.set_xlabel("# Reads", fontsize=15)
    ax.set_ylabel("Norm Difference per Peptide", fontsize=15)
    ax.hlines([0], 1, 50, linestyle=":", color="b")

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    fig.savefig("%s_normDiffBoxplots.png" % (opts.out), dpi=300, bbox_inches='tight')

    # Enriched peptide plot
    if opts.enrDir:

        #Read in enriched peptides from full dataset
        fullF = glob.glob("%s/Full~Full*txt" % (opts.enrDir))[0]
        fullEnrD = io.fileEmptyDict(fullF, header=False)

        allFs = glob.glob("%s/*.txt" % (opts.enrDir))

        x=[]
        total=[]
        true=[]
        false=[]

        for eachF in allFs:
            if eachF != fullF:
                thisEnrL = io.fileList(eachF, header=False)
                thisTrue = 0
                thisFalse = 0

                for p in thisEnrL:
                    if p in fullEnrD:
                        thisTrue+=1
                    else:
                        thisFalse+=1
        
                x.append(int(eachF.split("/")[-1].split("~")[0].split("_")[0]))
                total.append(len(thisEnrL))
                true.append(thisTrue)
                false.append(thisFalse)
    
        fig,ax = plt.subplots(figsize=(6,5),facecolor='w')

        ax.scatter(x, total, alpha=0.7, zorder=2, s=15, label="Total", color="b")
        ax.scatter(x, true, alpha=0.7, zorder=2, s=15, label="True", color="g")
        ax.scatter(x, false, alpha=0.7, zorder=2, s=15, label="False", color="r")
        ax.set_xlabel("# Reads", fontsize=15)
        ax.set_ylabel("Enriched peptides", fontsize=15)

        for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)

        ax.legend()

        fig.savefig("%s_enrichedPeps.png" % (opts.out), dpi=300, bbox_inches='tight')

#----------------------End of main()

def cumulativeNormDiff(fullD, subD):
    diff = 0
    for p, v in fullD.items():
        diff+=abs(v-subD[p])
    return diff

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

