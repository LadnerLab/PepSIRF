#! /usr/bin/env python3

# My custom libraries
import inout as io
import fastatools as ft
import kmertools as kt

# Built-in modules
import itertools as it
import glob, os, random, time, argparse
from collections import defaultdict

# 3rd party modules
import pandas as pd
import numpy as np
from scipy.signal import find_peaks

# This script is designed to read in lists of enriched peptides from high density PepSeq libraries (like RPTR1 and RPTR2) and output putative antibody epitopes

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-m', '--libMF',  help='A metadata file with info about the library peptides.', required=True)
	reqArgs.add_argument('-f', '--libFF', help='Fasta file with amino acid sequences for all peptides within the target library. ')
	reqArgs.add_argument('-p', '--plF', help='Tab-delimited file with info about the target sequences. Expected columns: "Source", "Protein", "Sequence"')
	reqArgs.add_argument('-o', '--out', help='Name for the ab-delimited output file that will be generated with information about the putative epitopes.')

	parser.add_argument('--minEpiSize', default=5, type=int, help='Minimum size of epitope to consider.')
	parser.add_argument('--maxEpiSize', default=30, type=int, help='Maximum size of epitope to consider. This should generally be the same as the size of the peptides. Otherwise, the epitope selected for individual enriched peptides will be somewhat arbitrary.')
	parser.add_argument('--maxEpisPerReg', default=2, type=int, help='Maximum number of epitopes per region to attempt to identify using the exhaustive approach.')
	parser.add_argument('--run2ndRound', default=False, action="store_true", help='Use this flag if you want to read in the inidivdual sample epitope calls and further process to the best dataset levels epitopes.')


	args = parser.parse_args()
	
	# Read in Library metadata
	libMD = pd.read_csv(args.libMF, sep="\t", header=0, index_col=0, keep_default_na=False)
	# Read in library peptide sequneces
	libFD = ft.read_fasta_dict(args.libFF)
	
	#Read in target sequence info
	plDF = pd.read_csv(args.plF, sep="\t", header=0)
	# Extract needed info about the target sequences
	plD, kpD = parseTargetSeqs(plDF, args.minEpiSize, args.maxEpiSize)

	# Read in enriched peptides and find epitopes
	epDir = "/Users/jtl276/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/DTRA_Projects/RAPTER/RPTR2/IM0206/aps/flexible/z10/aps-tsv/10Z-HDI95_0CS_75000raw"
	epFL = glob.glob(f"{epDir}/*_enriched.txt")
	
	allEpsD = {k:[] for k in ["Sample", "Virus", "Protein", "EpitopeSequence", "StartPos", "EndPos", "EpitopePositions", "Approach"]}
	
	# Step through the list of enriched peptides for each sample
	for fp in epFL:
		sampT = tuple(sorted(os.path.basename(fp)[:-13].split("~")))
		epL = io.fileList(fp, header=False)
		
		# Skip any samples that don't have enriched peptides
		if epL != [" "]:

			# Make sure that CtoS peptides are removed IF the WT version is also enriched
			# No double counting!
			# And, for the sake of simplicity, I'm also replacing CtoS with WT, if only CtoS is present
			c2s = [pn for pn in epL if libMD["CtoS"][pn]]
			toSwap = [pn for pn in c2s if libMD["CtoS"][pn] not in epL]
			for pn in toSwap:
				epL.append(libMD["CtoS"][pn])
			epL = [pn for pn in epL if pn not in c2s]
	
			# Generate starting lists (contained in a dictionary) for each protein with 0's for every position
			# And another dictionary to store info about which peptides overlap each position
			epD = {}
			pepsByPosD = {}
			# Stepping through each protein in the library
			for prot, l in plD.items():
				epD[prot] = [0]*l
				pepsByPosD[prot] = []
				for i in range(l):
					pepsByPosD[prot].append({})
	
			# Add info for each enriched peptide, at all overlapping positions
			for pn in epL:
				vir = libMD["Virus"][pn]
				prot = libMD["Protein"][pn]
				if libMD["Category"][pn] == "Focal":
					for i in range(int(libMD["StartPos"][pn]), int(libMD["EndPos"][pn])):
						epD[(vir,prot)][i]+=1
						pepsByPosD[(vir,prot)][i][pn]=""
			
			# Step through each library protein
			for protT, covL in epD.items():
				# Check to see if there are reactive peptides
				if sum(covL) > 0:
					failedRegs = findEpis(covL, sampT, allEpsD, protT, pepsByPosD, libMD, libFD, kpD, args)
					individual=0
					while len(failedRegs)>0:
						individual+=1
						newCovL = list(sum([np.array(x) for x in failedRegs]))
						failedRegs = findEpis(newCovL, sampT, allEpsD, protT, pepsByPosD, libMD, libFD, kpD, args, extraNote=f" after Individualx{individual}")
			
	# Write out putative epitopes
	epsDF = pd.DataFrame(allEpsD)
	epsDF.to_csv(args.out, sep="\t", index=False)
	
	if args.run2ndRound:
		comboEpsD = {k:[] for k in ["Sample", "Virus", "Protein", "EpitopeSequence", "StartPos", "EndPos", "EpitopePositions", "Approach"]}

		epFD = {}
		epD = {}
		pepsByPosD = {}
		# Stepping through each protein in the library
		for prot, l in plD.items():
			epD[prot] = [0]*l
			pepsByPosD[prot] = []
			for i in range(l):
				pepsByPosD[prot].append({})

		# Add info for each selected epitope, at all overlapping positions
		for j, row in epsDF.iterrows():
			epFD[j]=row["EpitopeSequence"]
			vir = row["Virus"]
			prot = row["Protein"]
			for i in row["EpitopePositions"].split("~"):
				epD[(vir,prot)][int(i)]+=1
				pepsByPosD[(vir,prot)][int(i)][j]=""

		# Step through each library protein
		for protT, covL in epD.items():
			# Check to see if there are reactive peptides
			if sum(covL) > 0:
				failedRegs = findEpis(covL, {"Combo"}, comboEpsD, protT, pepsByPosD, epsDF, epFD, kpD, args)
				individual=0
				while len(failedRegs)>0:
					individual+=1
					newCovL = list(sum([np.array(x) for x in failedRegs]))
					failedRegs = findEpis(newCovL, {"Combo"}, comboEpsD, protT, pepsByPosD, epsDF, epFD, kpD, args, extraNote=f" after Individualx{individual}")
					
		
		# Add column with information about which peptides fully contain the selected epitopes
		kmers2Peps = {}
		for n,s in libFD.items():
			if libMD["Category"][n]=="Focal":
				vir = libMD["Virus"][n]
				prot = libMD["Protein"][n]
				if (vir,prot) not in kmers2Peps:
					kmers2Peps[(vir,prot)] = defaultdict(list)
				
				seqKmers = []
				for k in range(args.minEpiSize, args.maxEpiSize+1):
					these = kt.kmerList(s,k)
					seqKmers+=these
					
					# Add in the WT kmers for CtoS peptides
					if libMD["CtoS"][n]:
						these = kt.kmerList(libFD[libMD["CtoS"][n]],k)
						seqKmers+=these
						
				for km in set(seqKmers):
					kmers2Peps[(vir,prot)][km].append(n)

		comboEpsD["Peptides"] = []
		for i,epseq in enumerate(comboEpsD["EpitopeSequence"]):
			theseEpPepsL = kmers2Peps[(comboEpsD["Virus"][i], comboEpsD["Protein"][i])][epseq]
			comboEpsD["Peptides"].append("~".join(theseEpPepsL))
		
		# Write out putative epitopes
		comboEpsDF = pd.DataFrame(comboEpsD)
		comboEpsDF = comboEpsDF.sort_values(by=['Virus', 'Protein', 'StartPos', 'EndPos'])
		comboEpsDF.to_csv(f"{args.out}_combo.tsv", sep="\t", index=False)
		
		# Combine overlapping epitopes 
		mergedRows = []
		rowCount=-1
		for i, row in comboEpsDF.iterrows():
			rowCount+=1
			if rowCount==0:
				priorRow = row
				priorPeps = set(row["Peptides"].split("~"))
				priorPos = set([int(j) for j in row["EpitopePositions"].split("~")])
				
			else:
				thesePeps = set(row["Peptides"].split("~"))
				thesePos = set([int(j) for j in row["EpitopePositions"].split("~")])
				if len(thesePos.intersection(priorPos))>0:
					priorRow["EpitopeSequence"] = f"{priorRow['EpitopeSequence']},{row['EpitopeSequence']}"
					priorRow["EndPos"] = row['EndPos']
					priorRow["EpitopePositions"] = "~".join([str(j) for j in list(range(priorRow["StartPos"],priorRow["EndPos"]))])
					priorRow["Approach"] = "Merged"
					priorPeps = thesePeps.union(priorPeps)
					priorPos = thesePos.union(priorPos)
					priorRow["Peptides"] = "~".join(sorted(list(priorPeps)))
					
				else:
					mergedRows.append(priorRow)
					priorRow = row
					priorPeps = thesePeps
					priorPos = thesePos
		mergedRows.append(priorRow)
		mergedDF = pd.DataFrame(mergedRows)
		mergedDF.to_csv(f"{args.out}_comboMerged.tsv", sep="\t", index=False)

#####-------------------->>>

def splitReactiveRegions(covL):
	regionsCL = regionCoords(covL)

	allSplit = {}
	for s,e in regionsCL:
		this = [0]*s + covL[s:e+1] + [0]*(len(covL)-(e+1))
		allSplit[(s,e)] = this

	return allSplit
		
def regionCoords(covL):
	regionsCL=[]

	inregion=0
	s=""
	e=""

	for i,c in enumerate(covL):
		if c > 0:
			if inregion:
				e=i
			else:
				inregion=1
				s=i
				e=i
		else:
			if inregion:
				regionsCL.append((s,e))
				inregion=0
				s=""
				e=""
	if inregion:
		regionsCL.append((s,e))
		
	return regionsCL

def parseTargetSeqs(plDF, minEpiSize, maxEpiSize):
	# Dictionary that will link a sequence to its length
	plD = {}
	# Dictionary that will link a sequences kmers to their positions in the sequence
	kpD={}
	
	for i,row in plDF.iterrows():
		plD[(row["Source"],row["Protein"])] = len(row["Sequence"])
		kpD[(row["Source"],row["Protein"])] = {}
		for ks in range(minEpiSize,maxEpiSize+1):
			for kmer,posL in kt.kmerDictPos(row["Sequence"],ks).items():
				kpD[(row["Source"],row["Protein"])][kmer]=posL
				
	return plD, kpD

def findEpis(covL, sampT, allEpsD, protT, pepsByPosD, libMD, libFD, kpD, args, extraNote=""):
	# To start, generate one coverage list that's specific to each reactive region
	# Positions overlapping other reactive regions will be converted to zeros
	regionsD = splitReactiveRegions(covL)
	
	# Variable to keep track of whether the epitope search has failed for some regions
	failed=[]
	
	# Step through each region and identify epitopes
	regNum=0
	for regCoords, regCov in regionsD.items():
		regNum+=1
		regPeps = []
		for i in range(regCoords[0], regCoords[1]+1):
			regPeps +=pepsByPosD[protT][i]
		regPeps = sorted(list(set(regPeps)))
		regSeqs = {pn:libFD[pn] for pn in regPeps}
		
		# Extract all relevant kmers from the peptides
		regKmers = set()
		regKmersBySeq = {}
		regKmers2Seq = defaultdict(list)
		regKmerScores = {}
		for n,s in regSeqs.items():
			seqKmers = []
			for k in range(args.minEpiSize, args.maxEpiSize+1):
				these = kt.kmerList(s,k)
				seqKmers+=these
			seqKmers = set(seqKmers)
			regKmersBySeq[s] = seqKmers
			regKmers.update(seqKmers)
			for km in seqKmers:
				regKmerScores[km] = sum([regCov[i] for i in kpD[protT][km]])
				regKmers2Seq[km].append(n)

		match=[]
		numEpi=0
		start = time.time()
		while len(match)==0 and numEpi<args.maxEpisPerReg:
			numEpi+=1
			for epis in it.combinations(regKmers, numEpi):
				
				episS = set(epis)
		
				# Set default to True, this will be turned to False once a seq without an overlap is encountered
				matchFound=True
				for n,s in regSeqs.items():
					thisOvlp = len(episS.intersection(regKmersBySeq[s]))
					if thisOvlp==0:
						matchFound=False
						break
					
				if matchFound:
					match.append(epis)
			# Check timing
			end = time.time()
		
		if match:
			# Calculate scores for each of the matching epi combos
			scoreD = defaultdict(list)
			for epis in match:
				epi_score = sum([regKmerScores[km] for km in epis])
				scoreD[epi_score].append(epis)
		
			# Select a set of epis with the highest score
			maxScore = max(scoreD.keys())
			chosen = random.choice(scoreD[maxScore])
			for ee in chosen:
				allEpsD["Sample"].append(",".join(sampT))
				allEpsD["Virus"].append(protT[0])
				allEpsD["Protein"].append(protT[1])
				allEpsD["EpitopeSequence"].append(ee)
				allEpsD["EpitopePositions"].append("~".join([str(x) for x in kpD[protT][ee]]))
				allEpsD["StartPos"].append(min(kpD[protT][ee]))
				allEpsD["EndPos"].append(max(kpD[protT][ee])+1)
				allEpsD["Approach"].append(f"Exhaustive{extraNote}")
		else:
			# Select the best epitope to start with by choosing the one that covers the greatest number of peptides, with the highest score
			numpepD = defaultdict(list)
			for km,pnL in regKmers2Seq.items():
				numpepD[len(pnL)].append(km)
			maxNP = max(numpepD.keys())
			
			scoreD = defaultdict(list)
			for km in numpepD[maxNP]:
				sc = regKmerScores[km]
				scoreD[sc].append(km)
			maxScore = max(scoreD.keys())
			
			# Add epitope to the dictionary
			chosen = random.choice(scoreD[maxScore])
			allEpsD["Sample"].append(",".join(sampT))
			allEpsD["Virus"].append(protT[0])
			allEpsD["Protein"].append(protT[1])
			allEpsD["EpitopeSequence"].append(chosen)
			allEpsD["EpitopePositions"].append("~".join([str(x) for x in kpD[protT][chosen]]))
			allEpsD["StartPos"].append(min(kpD[protT][chosen]))
			allEpsD["EndPos"].append(max(kpD[protT][chosen])+1)
			allEpsD["Approach"].append(f"Individual{extraNote}")
			
			# Modify covL and pepsByPosD to remove the coverage from peptides covered by the selected epitopes
			for pn in regKmers2Seq[chosen]:
				for i in range(int(libMD["StartPos"][pn]), int(libMD["EndPos"][pn])):
					regCov[i]-=1
					del(pepsByPosD[protT][i][pn])
					
#			print("Selected individual epitope based on coverage due to complex reactive region requiring more than the max specified epitopes:", chosen)

			failed.append(regCov)
			
	return failed

if __name__ == '__main__':
	main()
