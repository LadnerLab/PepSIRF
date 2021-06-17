#!/usr/bin/env python

import argparse, os, subprocess, glob
import inout as io
from collections import defaultdict    

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")
    parser.add_argument('--pepHeader', default="CodeName", help="Header name for metadata column with peptide names.")
    parser.add_argument('--tidHeader', default="SpeciesID", help="Header name for metadata column with species IDs.")
    parser.add_argument('--enrExt', default="_enriched.txt", help='Common file ending for files containing enriched sets of peptides. Once removed, the remaining filename should consist only of sample name(s)')
    parser.add_argument('-p', '--esPath', default = "~/Documents/GitHub/PepSIRF/extensions/enrichedScatterplots.py", help='Path to enrichedScatterplots.py.')
    parser.add_argument('--negMatrix',  help='Optional way to provide a separate data matrix for negative controls. If not provided, negative controls will be assumed to be from the same matrix as the experiemntal data.')
    parser.add_argument('-c', '--negControls',  help='Comma-delimited names of samples to use as negative controls. If this is not provided, then all of the samples in the negative control matrix will be used.')
    parser.add_argument('-i', '--negative_id', help='Optional approach for identifying negative controls. Provide a unique string at the start of all negative control samples.')

    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-e', '--enrDir',  help='Directory containing lists of enriched peptides.', required=True)
    reqArgs.add_argument('-t', '--tidFile',  help='File with list of taxIDs to explore.', required=True)
    reqArgs.add_argument("-d", "--data", help="Col-sum normalized read count file to pull data for plots.", required=True)
    reqArgs.add_argument('-m', '--meta', help='Metadata file that will be used to link peptide names to tax IDs.', required=True)
    reqArgs.add_argument('-o', '--outDir', help='Directory name for output files. Will be created.', required=True)

    args = parser.parse_args()

    #Create output dir
    genDir(args.outDir)

    # Read in tax IDs of interest
    tids = io.fileEmptyDict(args.tidFile, header=False)
    
    # Read in metadata
    pn2tid = io.fileDictHeader(args.meta, args.pepHeader, args.tidHeader)
    
    # Read in lists of enriched peptides
    enrFiles = glob.glob("%s/*%s" % (args.enrDir, args.enrExt))
    
    # Generate a dictionary to hold info for peptides of interest
    infD = {tid:defaultdict(list) for tid in tids}
    
    #Step through each enriched peptide file
    for eF in enrFiles:

        # Read in enriched peptides
        enrichedD = io.fileEmptyDict(eF, header=False)
        
        # Step through each enriched peptide and add to dictionary if it was designed from a taxon of interest
        for pn in enrichedD:
            thisTid = pn2tid[pn]
            if thisTid in tids:
                infD[thisTid][eF].append(pn)
    
    # Step through each tid of interest and generate output files
    for tid, hitInfo in infD.items():
        if len(hitInfo) > 0:
            
            tidEnrDir = "%s/%s" % (args.outDir, tid)
            
            # Create taxon-specific output dir
            genDir(tidEnrDir)
            
            enrFile = "%s/%s/%s_toPlot.tsv" % (args.outDir, tid, tid)
            
            # Write out sample-specific enriched peptide files
            for eF, eps in hitInfo.items():
                newPF = "%s/%s/%s" % (args.outDir, tid, os.path.basename(eF))
                io.writeList(eps, newPF)
            
            # Run enrichedScatterplots.py
            cmd = "%s -d %s -e %s --enrExt %s -o %s --plotLog 1 " % (args.esPath, args.data, tidEnrDir, args.enrExt, tidEnrDir)
            if args.negMatrix:
                 cmd += " --negMatrix %s " % (args.negMatrix)
            if args.negControls:
                 cmd += " --negControls %s " % (args.negControls)
            elif args.negative_id:
                 cmd += " --negative_id %s " % (args.negative_id)
            
            print(cmd)
            subprocess.run(cmd, shell=True)


#----------------------End of main()

def genDir(dirName):
    if os.path.isdir(dirName):
        print("The output diretory %s already exists!" % (dirName))
    else:
        os.mkdir(dirName)

    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

