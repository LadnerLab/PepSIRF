#!/usr/bin/env python

import argparse, os
import inout as io
#import matrixtools as mt
from collections import defaultdict
import numpy as np

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.interpolate import griddata
    
# Used to generate plots for identifying proper enrichment thresholds

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-x", "--xMat", help="Optionally, a separate score matricx can be provided for negative controls. If this option isn't used, negative control scores will be obtained from '--yMat'.")
    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")
    parser.add_argument('-o', '--outDir', default="thresholdScatterplots", help='Directory name for output files. Will be created.')
    parser.add_argument('--plotType', default="png", help='Type of plot to generate. This should be a file extension recognized by matplotlib, e.g., pdf, png, tiff.')

    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument("-y", "--yMat", help="Score matrix from which y-axis values will be derived.", required=True)
    reqArgs.add_argument("-s", "--samples", help="File containing sample names. One row should be provided per plot. Multiple names can be provided per line, separated by tabs; values from these samples will be average for the plot.", required=True)
    reqArgs.add_argument("-c", "--controls", help="Comma-separated list of control names. Average values from these samples will be plotted on the x-axis", required=True)
    reqArgs.add_argument("-t", "--thresh", help="Threshold file for plot contours. Tab-delimited file with one line per score matrix. The 1st column should contain a short descriptor that will be added to output names. The 2nd column should be the file path to the score matrix. And the 3rd column should be a comma-sep list of thresholds to contour on plots.")

    args = parser.parse_args()

    #Create output dir
    if os.path.isdir(args.outDir):
        print("The output diretory %s already exists!" % (args.outDir))
    else:
        os.mkdir(args.outDir)

    #Read in contour threshold file
    threshD = {}
    with open(args.thresh, "r") as fin:
        for line in fin:
            sampStr, filePath, ts = line.rstrip("\n").split("\t")
            threshD[sampStr] = [filePath, ts.split(",")]
    
    #Read in sample names
    sampleL = []
    with open(args.samples, "r") as fin:
        for line in fin:
            sampleL.append(line.rstrip("\n").split("\t"))
    
    # Read in score matrix for y-axis
    dataD = io.fileDictFullRowNames(args.yMat, args.delim)
    
    #Prep negative control data for plotting on x-axis, this will be the same across all of the plots
    if args.xMat:
        negD = io.fileDictFullRowNames(args.xMat, args.delim)
    else:
        negD = dataD

    negNames = args.controls.split(",")
    
    peptideNames = list(negD[negNames[0]].keys())
    
    x = np.array([np.mean([float(negD[sn][pn]) for sn in negNames]) for pn in peptideNames])
    xlog = np.log10(x+1)
    
    #Step through each score matrix to be used for drawing contours
    for scoreName, inf in threshD.items():
        scoreD = io.fileDictFullRowNames(inf[0], args.delim)
        zThresh = inf[1]

        #Step through each sample to plot
        for sNames in sampleL:

            #Average values for y-axis
            y = np.array([np.mean([float(dataD[sn][pn]) for sn in sNames]) for pn in peptideNames])
            ylog = np.log10(y+1)
            
            # Base heatmap
            heatmap, xedges, yedges = np.histogram2d(xlog, ylog, bins=(50,70))
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            heatmap.T[heatmap.T == 0.0] = np.nan

            # Z values for contours
            z = np.array([np.mean([float(scoreD[sn][pn]) for sn in sNames]) for pn in peptideNames])
            points = np.column_stack((np.log10(x+1), np.log10(y+1)))

            grid_x, grid_y = np.mgrid[0:max(xlog):0.01, 0:max(ylog):0.01]

            grid_z = griddata(points, z, (grid_x, grid_y), method='linear')

            #Make plot
            fig,ax = plt.subplots(1,1,figsize=(7,7),facecolor='w')
            ax.imshow(heatmap.T, extent=extent, origin='lower', cmap="viridis_r")

            CS = ax.contour(grid_x, grid_y, grid_z, levels=zThresh, colors='k')
            ax.clabel(CS, inline=1, fontsize=14, fmt='%d')

            ax.set_xlim(xedges[0]-0.1, xedges[-1]+0.1)
            ax.set_ylim(yedges[0]-0.1, yedges[-1]+0.1)
            
            fig.savefig("%s/%s_%s.%s" % (args.outDir, "~".join(sNames), scoreName, args.plotType), dpi=300, bbox_inches='tight')
            plt.close(fig)
            
#----------------------End of main()

def makeDirName(args):
    dirName = ""
    infD = io.fileDict(args.thresh, header=False)
    for k, v in infD.items():
        typ = k.split("_")[-1].split(".")[-2]
        val = v.replace(",", "-")
        dirName+="%s%s_" % (val, typ)
    if args.raw:
        dirName+="%sraw" % (args.rawThresh.replace(",", "-"))
        return dirName
    else:
        return dirName[:-1]

def boxplot(counts, outname, thresh = None):
    fig,ax = plt.subplots(1,1,figsize=(4, 4),facecolor='w')
    ax.boxplot([float(x) for x in counts])
    if thresh:
        ax.hlines([float(x) for x in thresh.split(",")], 0.5, 1.5, linestyle="--")
    plt.savefig(outname,dpi=300,bbox_inches='tight')
    
    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

