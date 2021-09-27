#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import argparse, glob, os
import pandas as pd
import seaborn as sns
import numpy as np
import inout as io
import matrixtools as mt
#import statistics as stat
#import itertools as it
#from collections import defaultdict



# This script reads in 1) Z scores and 2) pep map files produced by the PepSIRF deconv module
# and produces boxplots of the Z score distributions for each species called as seropositive

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-z', '--zMat',  help='Z score matrix for generating boxplots.', required=True)
    reqArgs.add_argument('-m', '--mapDir', help='Name of directory containing pep maps.', required=True)
    reqArgs.add_argument('-e', '--mapExt',  help='Common file extension for all peptide map files.', required=True)
    reqArgs.add_argument('-o', '--outDir',  help="Directory for output files. It will be created if it doesn't already exist", required=True)

    parser.add_argument('--meta',  help='Tab-delimited metadata file with columns linking species ID to species name')
    parser.add_argument('--sid', default="SpeciesID", help='Header in metadata file for column with species ID.')
    parser.add_argument('--sp',  default="Species", help='Header in metadata file for column with species name.')
    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")

    parser.add_argument('-x', '--xHead', default="SpeciesID", help='Header in generated data files for x-axis boxplot categories.')
    parser.add_argument('-y', '--yHead', default="Zscore", help='Header in generated data files for y-axis of boxplots.')

    parser.add_argument('--hue',  help='Header in data file for "hue", or subcategories for x-axis.')
    parser.add_argument('--xLab', help="String for x-label. If not provided, --xHead is used.")
    parser.add_argument('--yLab', help="String for y-label. If not provided, --yHead is used.")
    parser.add_argument('--yLim', help="Comma-delim floats to use as the min and max y-axis limits.")
    parser.add_argument('--xRotate', default=False, action="store_true", help="Use this option if you want the x-labels to be rotated vertically.")
    parser.add_argument('--noStripPlot', default=False, action="store_true", help="Use this option if you want to only plot the boxes, not the corresponding strip plots.")
    parser.add_argument('--width', default=5, type=int, help="Figure width.")
    parser.add_argument('--height', default=4, type=int, help="Figure height.")
    parser.add_argument('--logY', default=False, action="store_true", help="Use this option if you want the y-axis to be coversted to log scale. If some of the values are <=0, an integer will be added to all y-axis values prior to conversion")
#     parser.add_argument('--include', help="Header,Value pairs used to indicate a subset of rows to include.", nargs='*')
#     parser.add_argument('--exclude', help="Header,Value pairs used to indicate a subset of rows to exclude.", nargs='*')
    parser.add_argument('--axhline', help="Y-axis value at which a horizontal line should be drawn.")

    opts = parser.parse_args()

    #If provided, create linkage map between species ID and name. 
    if opts.meta:
        sid2sp = io.fileDictHeader(opts.meta, opts.sid, opts.sp)
        print("Read in species ID/name linkage file.")
    else:
        sid2sp = {}
    
    # Read in Z score data file
    zD = mt.parseCounts(opts.zMat, delim=opts.delim)
    print("Read in Zscore data matrix.")

    # Generate output directory, if it doesn't already ecist
    if not os.path.isdir(opts.outDir):
        os.mkdir(opts.outDir)
    else:
        print("Warning: %s already exists" % (opts.outDir))
    
    # Step through each peptide map file
    mapFiles = glob.glob("%s/*%s" % (opts.mapDir, opts.mapExt))
    print("Processing %d map files." % (len(mapFiles)))
    
    boxInfoFiles=[]
    for mf in mapFiles:
        samples = mf.split("/")[-1][:-len(opts.mapExt)].split("~")
        mapInfoD = io.fileDictHeader(mf, "Peptide", "Assigned Ids")
        if sum([1 for s in samples if s in zD]) == len(samples):
            thisOut = "%s/%s_zScores.tsv" % (opts.outDir, "~".join(samples))
            boxInfoFiles.append(thisOut)
            with open(thisOut, "w") as fout:
                fout.write("SpeciesID\tSpecies\tReplicate\tCodeName\tZscore\n")
                for pid, sid in mapInfoD.items():
                    for s in samples:
                        
                        try:
                            sp = sid2sp[sid]
                        except:
                            sp = ""
                        
                        fout.write("%s\t%s\t%s\t%s\t%.3f\n" % (sid, sp, s, pid, zD[s][pid]))
        else:
            print("One or more of these samples were not found in the Z score matrix: %s" % ("\t".join(samples)))

    #Clear out memory by deleting the large Z score matrix after done using
    del(zD)

    if len(boxInfoFiles) > 0:
        for bif in boxInfoFiles:
            dataD = io.fileDictFull(bif, delim=opts.delim)
            dataD[opts.yHead] = [float(a) for a in dataD[opts.yHead]]

            
            if opts.logY:
                opts.addInt=False
                minY = min(dataD[opts.yHead])
                if minY<=0:
                    opts.addInt = 1
                    while opts.addInt + minY <=0:
                        opts.addInt+=1
                    dataD[opts.yHead] = [a+opts.addInt for a in dataD[opts.yHead]]
                dataD[opts.yHead] = [np.log10(a) for a in dataD[opts.yHead]]
    
            df = pd.DataFrame(dataD)
    
        #     if opts.include:
        #         for each in opts.include:
        #             header, val = each.split(",")
        #             df = df[df[header] == val]
        # 
        #     if opts.exclude:
        #         for each in opts.exclude:
        #             header, val = each.split(",")
        #             df = df[df[header] != val]
    
            catBoxplot(df, opts, colorHead=opts.hue, xLab=opts.xLab, yLab=opts.yLab, out="%s_boxes.png" % (bif[:-4]))

#----------------------End of main()

def catBoxplot(df, opts, colorHead=None, xLab=None, yLab=None, out=None):
    if xLab==None:
        xLab=opts.xHead
    if yLab==None:
        yLab=opts.yHead
        
    fig,ax = plt.subplots(1,1,figsize=(opts.width, opts.height),facecolor='w')
    
    if colorHead:
        sns.boxplot(x=opts.xHead, y=opts.yHead, hue=colorHead, data=df, ax=ax, fliersize=0, zorder=2)
        if not opts.noStripPlot:
            sns.stripplot(x=opts.xHead, y=opts.yHead, hue=colorHead, data=df, jitter=True, dodge=True, linewidth=0.5, ax=ax, zorder=3)
    
    else:
        sns.boxplot(x=opts.xHead, y=opts.yHead, data=df, ax=ax, fliersize=0, zorder=2)
        if not opts.noStripPlot:
            sns.stripplot(x=opts.xHead, y=opts.yHead, data=df, jitter=True, dodge=True, linewidth=0.5, ax=ax, zorder=3)
    
    if opts.axhline:
        ax.axhline(float(opts.axhline), 0, 1, zorder=1, ls='dotted', c='k')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax.set_xlabel(xLab, fontsize=10)
    if opts.logY:
        if opts.addInt:
            ax.set_ylabel("log10(%s+%d)" % (yLab, opts.addInt), fontsize=20)
        else:
            ax.set_ylabel("log10(%s)" % yLab, fontsize=20)
    else:
        ax.set_ylabel(yLab, fontsize=20)
    ax.tick_params(labelsize=15)
    
    if opts.xRotate:
        plt.xticks(rotation="vertical")
    
    if opts.yLim:
        yMin, yMax = [float(y) for y in opts.yLim.split(",")]
        ax.set_ylim(yMin, yMax)
    
    if out:
        plt.savefig(out,dpi=300,bbox_inches='tight')
    
    plt.close()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

