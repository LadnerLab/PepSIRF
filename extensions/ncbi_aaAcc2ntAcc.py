#!/usr/bin/env python3

import entrezpy.conduit

### WARNING: To be able to access the result from entrezpy, I had to do the following:
### 1. I manually opened this file: /Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages/entrezpy/efetch/efetch_analyzer.py
### 2. I added this new line of code following the print statement on Line 56: self.output = self.norm_response(response, request.rettype).rstrip('\n')
### 3. I commented out the print statement on Line 56
### 4. I saved the file with these changes

import argparse, os
import pandas as pd
# import fastatools as ft # Available from https://github.com/jtladner/Modules
# import inout as io		 # Available from https://github.com/jtladner/Modules
# from collections import defaultdict

#This script is designed to generate subset fasta files to serve as inputs for sequence-based clustering


def main():
	parser = argparse.ArgumentParser(description="Used to fetch NCBI nt accessions corresponding to associated NCBI aa accessions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	# Required arguments
	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument("-i", "--input", metavar="", required=True, help="Tab-delimited file with at least one column that includes NCBI amino acid accession numbers. Extra columns are fine and will be maintained, as is, in the output.")
	reqArgs.add_argument("-o", "--output", metavar="", required=True, help="Name for tab-delimited output file, which will be created from scratch or written over, if it already exists.")
	reqArgs.add_argument("-e", "--email", required=True, help="User's email address, in case NCBI needs to contact the user.")

	# Optional arguments
	parser.add_argument("--aaHeader", default="NCBIaccession-AA", help="Header within input file for the column that includes the amino acid accessions.")
	parser.add_argument("--ntHeader", default="NCBIaccession-NT", help="Header to be used in putput file and/or within the input file for the column that includes the nucleotide accessions.")
	parser.add_argument("--idDelim", default=",", help="Delimiter separating distinct amino acid accessions.")
	parser.add_argument("--colDelim", default="\t", help="Delimiter separating distinct columns for input and output files.")
	parser.add_argument("--tries", type=int, default=5, help="Number of numbers to attempt to retrieve an nt accession, for each aa accession.")
	parser.add_argument("--reportInc", type=int, default=1000, help="Print a message to the screen after processing this number of input lines.")

	args=parser.parse_args()

	#Read in input file
	df = pd.read_csv(args.input, sep=args.colDelim, header=0, dtype=str,keep_default_na=False)
	#Verify that expected header is present, only continue if it is
	if args.aaHeader not in df.columns:
		print(f"{args.aaHeader} is not a column name in {args.input}. Please revise your input file and/or command.")
	else:
		#Initiate a counter for tracking
		counter=0
		
		#Open the output file for writing
		with open(args.output, "w") as fout:
			# Write header for output
			headL = list(df.columns)
			
			# Check to see if the input file already has a column for the nt accessions, and if not, add it for the output
			if args.ntHeader in headL:
				inputHasNt=True
			else:
				inputHasNt=False
				headL.append(args.ntHeader)
			
			# Write out header
			headS = args.colDelim.join(headL)
			fout.write(f"{headS}\n")
			
			# Step through each row in the input
			for i,row in df.iterrows():
				aaAccL = row[args.aaHeader].split(args.idDelim)
				if inputHasNt:
					ntAccL = row[args.ntHeader].split(args.idDelim)
					for j,pid in enumerate(aaAccL):
						if ntAccL[j] == "":
							ntAccL[j] = prot2nucACC(pid, email=args.email, timesToTry = args.tries)
				else:
					ntAccL = [prot2nucACC(pid, email=args.email, timesToTry = args.tries) for pid in aaAccL]
					
				#Write out new line
				ntAccS = args.idDelim.join(ntAccL)
				outL = []
				for head in headL:
					if head==args.ntHeader:
						outL.append(ntAccS)
					else:
						outL.append(row[head])
				outS = args.colDelim.join(outL)
				fout.write(f"{outS}\n")
				
				#Increment counter
				counter+=1
				if counter%args.reportInc==0:
					print(f"{counter} lines processed.")
	
	
###---------------End of main()-----------------------------------

# Function for converting an AA accession into an NT accession
def prot2nucACC(pid, email, timesToTry = 1):
	if len(pid)>=2:
		# In case the version number is missing, add ".1"
		if pid[-2]!=".":
			pid+=".1"
	
		c = entrezpy.conduit.Conduit(email)   # Create new Conduit instance
		fetch_ppl = c.new_pipeline()		  # Create empty Conduit pipeline 
		sid = fetch_ppl.add_link({'dbfrom': 'protein', 'db': 'nuccore', 'id':pid})
		fid = fetch_ppl.add_fetch({'retmode': 'text', 'rettype': 'acc'}, dependency=sid)
		
		complete = False
		counter = 0
		
		while counter<timesToTry and complete==False:
			counter+=1
			try:
				c.run(fetch_ppl)
				complete = True
			except Exception:
				complete = False
		
		# Return the nuc accession
		if complete:
			try: 
				return c.analyzers[fid].output
			except Exception:
				return ""
		else:
			return ""
	else:
		return ""

if __name__ == "__main__":
	main()



