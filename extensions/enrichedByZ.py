#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import argparse, os

import inout as io
import pandas as pd
import numpy as np
from collections import defaultdict

import seaborn as sns

#colors = ["#D55E00", "#E69F00", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#999999"]
#colors = ["#e41a1c", "#377eb8", "#984ea3", "#4daf4a", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-d', '--dataF',  help='A tab-delimited file Z score data.', required=True)
	reqArgs.add_argument('-o', '--out', help='Basename for output files.', required=True)
	
	# Flags associated with running this either in batch mode or one protein at a time
	parser.add_argument('-m', '--smF',  help='A tab-delimited metadata file with info about the samples of interest. One of these columns should define categories of interest. The first column should contain a sample identifier.')
#	reqArgs.add_argument('-c', '--compCat', help='Sample metadata column to use for the categories that will be compared on a single plot.', required=True)

#	parser.add_argument('--colorMap', help='Optional tab-delimited file linking compCat categories to colors to use for the graphs. There should NOT be a header row. First column = catgroy name; second column = color hexcode or name.')
	parser.add_argument('--zRange', default="10,1000,1", help='Three comma-delim values indicating the range of Z scores to use for generating plots: min, max, stepsSize.')
#	parser.add_argument('--plotCat', help='Optional sample metadata column to use for splitting samples into multiple plots.')
	parser.add_argument('--logTrans', default=False, action="store_true", help='Use this flag if you want to log transform the Z scores for plotting and calculation of AUC.')
#	parser.add_argument('--noFigs', default=False, action="store_true", help='Use this flag if you do not want to generate the associated figures.')
	parser.add_argument('--plotWidth', default=6, type=int, help='Width of output boxplot, if generated.')
	parser.add_argument('--plotHeight', default=3.5, type=int, help='Height of output boxplot, if generated.')

	args = parser.parse_args()
	
	# Read in Z score data
	df = pd.read_csv(args.dataF, sep="\t", header=0, index_col=0)

	# Generate range of Z scores to use for the analysis
	args.zRange = [int(i) for i in args.zRange.split(",")]
	zRange = range(args.zRange[0], args.zRange[1], args.zRange[2])

	# Read sample metadata
#	smDF = pd.read_csv(args.smF, sep="\t", header=0, index_col=0, keep_default_na=False, dtype=str)

	# Generate data for covByZ statistic and format data for input to Seaborn for plotting
	longData = {"Sample":[], "NumReact":[], "zThresh":[],} #"Category":[], "PlotCategory":[]}
	
	for sn in df.columns:
		zL = list(df[sn])
		for zT in zRange:
			zL = [z for z in zL if z>=zT]
			longData["Sample"].append(sn)
			longData["NumReact"].append(len(zL))
#			longData["Category"].append(cat)
#			longData["PlotCategory"].append(plotCat)
			if args.logTrans:
				longData["zThresh"].append(np.log2(zT))
			else:
				longData["zThresh"].append(zT)
			
		
	longDF = pd.DataFrame(longData)
	
	areaD = {}
	
	# Calculate areas
	# Step through each individual sample
	snL = sorted(list(set(longDF["Sample"])))
	for sn in snL:
		sampDF = longDF[longDF["Sample"]==sn]
		sampDF = sampDF.sort_values(by=['zThresh'])
			
		# Calculate Area
		area = np.trapezoid(sampDF['NumReact'], sampDF['zThresh'])
		areaD[sn] = f"{area:.4f}"

	# Generate output file
	if args.logTrans:
		outscores = f"{args.out}_enrichedByZ_AUC-log2.tsv"
	else:
		outscores = f"{args.out}_enrichedByZ_AUC.tsv"

	with open(outscores, "w") as fout:
		fout.write(f"Sample\tScore\n")
		
		# Step through each sample
		for sn, strA in areaD.items():
			fout.write(f"{sn}\t{strA}\n")


# 		# Generate plot, if requested
# 		if not args.noFigs:
# 			fig,ax = plt.subplots(1,1,figsize=(args.plotWidth, args.plotHeight),facecolor='w')
# 			
# 			lp = sns.lineplot(data=plotDF, x="zThresh", y="PropReact", hue="Category", palette=colorD)    
# 			sns.move_legend(lp, "upper left", bbox_to_anchor=(1, 1))
# 			
# 			ax.set_ylabel('Proportion Reactive', fontsize=16)
# 			if args.logTrans:
# 				ax.set_xlabel('Minimum Z score (log2)', fontsize=16)
# 			else:
# 				ax.set_xlabel('Minimum Z score', fontsize=16)
# 		
# 				
# 			plt.grid(linestyle = '--', linewidth = 0.5)
# 	# 		ax.set_ylim(0, 0.6)
# 			outname = f"{args.out}_propReactRibbon_{tp}"
# 			if args.logTrans:
# 				outname+="_log2"
# 			fig.savefig(f"{outname}.pdf", bbox_inches="tight")
# 			fig.savefig(f"{outname}.png", bbox_inches="tight", dpi=300)


######------------------------------------------------------------------------>>>>>>>>>>>>


if __name__ == '__main__':
	main()

