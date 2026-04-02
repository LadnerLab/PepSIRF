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

colors = ["#D55E00", "#E69F00", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#999999"]
#colors = ["#e41a1c", "#377eb8", "#984ea3", "#4daf4a", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-m', '--smF',  help='A tab-delimited metadata file with info about the samples of interest. One of these columns should define categories of interest. The first column should contain a sample identifier.', required=True)
	reqArgs.add_argument('-d', '--dataF',  help='A tab-delimited file Z score data.', required=True)
	reqArgs.add_argument('-c', '--compCat', help='Sample metadata column to use for the categories that will be compared on a single plot.', required=True)
	reqArgs.add_argument('-o', '--out', help='Basename for output files.', required=True)
	
	# Flags associated with running this either in batch mode or one protein at a time
	parser.add_argument('-b', '--batchF',  help='Path to tab-delimited file with no header and two columns: protein name, associated map file in format described for --mapF. Using this flag is an alternative to providing --mapF. Each row in the input file will be equivalent to one command mannually specifying --mapF.', required=False)
	parser.add_argument('-p', '--mapF',  help='Path to file linking peptides to positions in target protein. This file should have two tab-delim columns: ProbeName, AlignPos. ProbeName should contain peptide names, while AlignPos should be ~ delimited positions.', required=False)

	parser.add_argument('--colorMap', help='Optional tab-delimited file linking compCat categories to colors to use for the graphs. There should NOT be a header row. First column = catgroy name; second column = color hexcode or name.')
	parser.add_argument('--zRange', default="10,500,1", help='Three comma-delim values indicating the range of Z scores to use for generating plots: min, max, stepsSize.')
	parser.add_argument('--plotCat', help='Optional sample metadata column to use for splitting samples into multiple plots.')
	parser.add_argument('--logTrans', default=False, action="store_true", help='Use this flag if you want to log transform the Z scores for plotting and calculation of AUC.')
	parser.add_argument('--noFigs', default=False, action="store_true", help='Use this flag if you do not want to generate the associated figures.')
	parser.add_argument('--plotWidth', default=6, type=int, help='Optional sample metadata column to use for splitting samples into multiple plots.')
	parser.add_argument('--plotHeight', default=3.5, type=int, help='Optional sample metadata column to use for splitting samples into multiple plots.')

	args = parser.parse_args()
	
	# Check to determine what mode to run the script in Batch or Single
	if args.batchF:
		protD = io.fileDict(args.batchF, header=False)
		print(f"Running in Batch Mode. Processing {len(protD)} proteins.")
	elif args.mapF:
		print(f"Running in Standard Mode.")
		if args.logTrans:
			protD = {"AUC_log2": args.mapF}
		else:
			protD = {"AUC": args.mapF}
	else:
		print("ERROR: You must either use '--batchF' or '--mapF' flag.")
		return False

	# Read in Z score data
	df = pd.read_csv(args.dataF, sep="\t", header=0, index_col=0)
	# Generate range of Z scores to use for the analysis
	args.zRange = [int(i) for i in args.zRange.split(",")]
	zRange = range(args.zRange[0], args.zRange[1], args.zRange[2])


	# Define colors to use for plot
	if not args.noFigs:
		if args.colorMap:
			colorD = io.fileDict(args.colorMap, header=False)
		else:
			colorD = {}
			for i, c in enumerate(sorted(list(set(longData["Category"])))):
				colorD[c] = colors[i]

	# Read sample metadata
	smDF = pd.read_csv(args.smF, sep="\t", header=0, index_col=0, keep_default_na=False, dtype=str)

	
	outStrD = {}
	protColumns = sorted(list(protD.keys()))

	pCount=0
	for outN in protColumns:
		mapF = protD[outN]
		# Read in alignment map
		try:
			mapD = io.fileDictHeader(mapF, "ProbeName", "AlignPos")
			mapD = {k:[int(p) for p in v.split("~")] for k,v in mapD.items()}
		except:
			print(f"Skipping {outN}. There was an issue reading in map file.")
			continue

		pCount+=1
		
		# Find number of unique positions in map
		allPos = []
		for v in mapD.values():
			allPos+=v
		numUniqPos = len(set(allPos))
	
		# Generate data for covByZ statistic and format data for input to Seaborn for plotting
		longData = {"Sample":[], "Category":[], "PropReact":[], "zThresh":[], "PlotCategory":[]}
		
		for sn in df.columns:
			if sn in smDF[args.compCat]:
				cat = smDF[args.compCat][sn]
				if args.plotCat:
					plotCat = smDF[args.plotCat][sn]
				else:
					plotCat = "Main"
				
				prevCov = -99999
				for zT in zRange:
					if prevCov==0:
						propCov=0
					else:
						posCov = posCovered(df[sn], zT, mapD)
						propCov = posCov/numUniqPos
						prevCov = posCov
						
					longData["Sample"].append(sn)
					longData["Category"].append(cat)
					longData["PropReact"].append(propCov)
					longData["PlotCategory"].append(plotCat)
					if args.logTrans:
						longData["zThresh"].append(np.log2(zT))
					else:
						longData["zThresh"].append(zT)
			
			elif pCount==1:
				print(f"Skipping {sn}, not in sample metadata")
				
			
		longDF = pd.DataFrame(longData)

		for tp in sorted(list(set(longData["PlotCategory"]))):

			plotDF = longDF[longDF["PlotCategory"]==tp]

			# Calculate areas
			# Step through each individual sample
			snL = sorted(list(set(plotDF["Sample"])))
			for sn in snL:
				if sn not in outStrD:
					if args.plotCat:
						outStrD[sn] = [f"{smDF[args.compCat][sn]}\t{tp}"]
					else:
						outStrD[sn] = [f"{smDF[args.compCat][sn]}"]
					
				sampDF = plotDF[plotDF["Sample"]==sn]
				sampDF = sampDF.sort_values(by=['zThresh'])
					
				# Calculate Area
				area = np.trapezoid(sampDF['PropReact'], sampDF['zThresh'])
				outStrD[sn].append(f"{area:.4f}")


			# Generate plot, if requested
			if not args.noFigs:
				fig,ax = plt.subplots(1,1,figsize=(args.plotWidth, args.plotHeight),facecolor='w')
				
				lp = sns.lineplot(data=plotDF, x="zThresh", y="PropReact", hue="Category", palette=colorD)    
				sns.move_legend(lp, "upper left", bbox_to_anchor=(1, 1))
				
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
		
		if pCount % 100 == 0:
			print(f"Processing {pCount/len(protColumns)*100}% complete.")
		
	# Generate output file
	if args.logTrans:
		outscores = f"{args.out}_covByZ_AUC-log2.tsv"
	else:
		outscores = f"{args.out}_covByZ_AUC.tsv"

	with open(outscores, "w") as fout:
		protColStr = "\t".join(protColumns)
		if args.plotCat:
			fout.write(f"Sample\tGroup\tCategory\t{protColStr}\n")
		else:
			fout.write(f"Sample\tGroup\t{protColStr}\n")
		
		# Step through each sample
		for sn, strL in outStrD.items():
			outS = "\t".join(strL)
			fout.write(f"{sn}\t{outS}\n")


######------------------------------------------------------------------------>>>>>>>>>>>>

def posCovered(dSeries, zT, mapD):
	posCov = []
	meetThresh = dSeries[dSeries>=zT]
	for pn, val in meetThresh.items():
		if pn in mapD:
			posCov+=mapD[pn]
	return len(set(posCov))


if __name__ == '__main__':
	main()

