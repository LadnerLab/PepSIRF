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

colors = ["#D55E00", "#E69F00", "CC79A7", "#009E73", "F0E442", "#0072B2", "E69F00", "#999999"]

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-p', '--mapF',  help='Path to file linking peptides to positions in target protein. This file should have two tab-delim columns: ProbeName, AlignPos. ProbeName should contain peptide names, while AlignPos should be ~ delimited positions.', required=True)
	reqArgs.add_argument('-m', '--smF',  help='A tab-delimited metadata file with info about the samples of interest. One of these columns should define categories of interest. The first column should contain a sample identifier.', required=True)
	reqArgs.add_argument('-d', '--dataF',  help='A tab-delimited file Z score data.', required=True)
	reqArgs.add_argument('-c', '--compCat', help='Sample metadata column to use for the categories that will be compared on a single plot.', required=True)
	reqArgs.add_argument('-o', '--out', help='Basename for output files.', required=True)

	parser.add_argument('--colorMap', help='Optional tab-delimited file linking compCat categories to colors to use for the graphs. There should NOT be a header row. First column = catgroy name; second column = color hexcode or name.')
	parser.add_argument('--zRange', default="10,500,1", help='Three comma-delim values indicating the range of Z scores to use for generating plots: min, max, stepsSize.')
	parser.add_argument('--plotCat', help='Optional sample metadata column to use for splitting samples into multiple plots.')
	parser.add_argument('--logTrans', default=False, action="store_true", help='Use this flag if you want to log transform the Z scores for plotting and calculation of AUC.')
	parser.add_argument('--plotWidth', default=6, type=int, help='Optional sample metadata column to use for splitting samples into multiple plots.')
	parser.add_argument('--plotHeight', default=3.5, type=int, help='Optional sample metadata column to use for splitting samples into multiple plots.')

	args = parser.parse_args()

	# Read in alignment map
	mapD = io.fileDictHeader(args.mapF, "ProbeName", "AlignPos")
	mapD = {k:[int(p) for p in v.split("~")] for k,v in mapD.items()}
	
	# Find number of unique positions in map
	allPos = []
	for v in mapD.values():
		allPos+=v
	numUniqPos = len(set(allPos))
	
	# Read sample metadata
	smF = '/Users/jtl276/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/DTRA_Projects/RAPTER/RPTR2/IM0217_IMO218/updated_separate_aps/JTL_analysis/sampleLevel_metadata.tsv'
	smDF = pd.read_csv(args.smF, sep="\t", header=0, index_col=0, keep_default_na=False, dtype=str)
	
	
	# Read in Z score data
	df = pd.read_csv(args.dataF, sep="\t", header=0, index_col=0)
	
	# Generate range of Z scores to use for the analysis
	args.zRange = [int(i) for i in args.zRange.split(",")]
	zRange = range(args.zRange[0], args.zRange[1], args.zRange[2])
	

	# Format data for input to Seaborn for plotting
	longData = {"Sample":[], "Category":[], "PropReact":[], "zThresh":[], "PlotCategory":[]}
	
	for sn in df.columns:
		if sn in smDF[args.compCat]:
			cat = smDF[args.compCat][sn]
			if args.plotCat:
				plotCat = smDF[args.plotCat][sn]
			else:
				plotCat = "Main"
			
			for zT in zRange:
				posCov = posCovered(df[sn], zT, mapD)
				propCov = len(posCov)/numUniqPos
				longData["Sample"].append(sn)
				longData["Category"].append(cat)
				longData["PropReact"].append(propCov)
				longData["PlotCategory"].append(plotCat)
				if args.logTrans:
					longData["zThresh"].append(np.log2(zT))
				else:
					longData["zThresh"].append(zT)
		
		else:
			print(f"Skipping {sn}, not in sample metadata")
		
	longDF = pd.DataFrame(longData)


	# Define colors to use for plot
	if args.colorMap:
		colorD = io.fileDict(args.colorMap, header=False)
	else:
		colorD = {}
		for i, c in enumerate(sorted(list(set(longData["Category"])))):
			colorD[c] = colors[i]

	
	if args.logTrans:
		outscores = f"{args.out}_covByZ_AUC-log2.tsv"
	else:
		outscores = f"{args.out}_covByZ_AUC.tsv"
	
	with open(outscores, "w") as fout:
		if args.logTrans:
			fout.write("Group\tCategory\tSample\tAUC_log2\n")
		else:
			fout.write("Group\tCategory\tSample\tAUC\n")

		# Generate plot
		for tp in sorted(list(set(longData["PlotCategory"]))):
		
			plotDF = longDF[longDF["PlotCategory"]==tp]
		
			fig,ax = plt.subplots(1,1,figsize=(args.plotWidth, args.plotHeight),facecolor='w')
			
			sns.lineplot(data=plotDF, x="zThresh", y="PropReact", hue="Category", palette=colorD)    
			
			ax.set_ylabel('Proportion Reactive', fontsize=16)
			if args.logTrans:
				ax.set_xlabel('Minimum Z score (log2)', fontsize=16)
			else:
				ax.set_xlabel('Minimum Z score', fontsize=16)
		
				
			plt.grid(linestyle = '--', linewidth = 0.5)
	# 		ax.set_ylim(0, 0.6)
			outname = f"{args.out}_propReactRibbon_{tp}"
			if args.logTrans:
				outname+="_log2"
			fig.savefig(f"{outname}.pdf", bbox_inches="tight")
			fig.savefig(f"{outname}.png", bbox_inches="tight", dpi=300)
			
			# Calculate areas
			# Step through each individual sample
			snL = sorted(list(set(plotDF["Sample"])))
			for sn in snL:
				sampDF = plotDF[plotDF["Sample"]==sn]
				sampDF = sampDF.sort_values(by=['zThresh'])
				
			# Calculate Area
				area = np.trapezoid(sampDF['PropReact'], sampDF['zThresh'])
				fout.write(f"{smDF[args.compCat][sn]}\t{tp}\t{sn}\t{area}\n")
	
######------------------------------------------------------------------------>>>>>>>>>>>>

def posCovered(dSeries, zT, mapD):
	posCov = []
	meetThresh = dSeries[dSeries>=zT]
	for pn, val in meetThresh.items():
		if pn in mapD:
			posCov+=mapD[pn]
	return sorted(list(set(posCov)))


if __name__ == '__main__':
	main()

