#!/usr/bin/env python

import argparse, os, glob
import inout as io
import pandas as pd
from scipy.stats import poisson
pd.options.mode.copy_on_write = True
from collections import defaultdict as dd
# import numpy as np
# import inout as io


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-e', '--enrDir',  help='Directory containing enriched peptide files.', required=True)
    reqArgs.add_argument('-r', '--protMeta',  help='Protein metadata file.', required=True)
    reqArgs.add_argument('-p', '--pepMeta',  help='Peptide metadata file.', required=True)
    reqArgs.add_argument('-o', '--outDir', help='Name/path for directory to place output files within. It will be created if it does not already exist', required=True)

    parser.add_argument('--protNameHeader', default="Sequence",  help='Header for protein name in protein metadata file.')
    parser.add_argument('--taxHeader', default="MultiProteinCluster",  help='Header for taxonomic category name in protein metadata file.')
    parser.add_argument('--pepNameHeader', default="CodeName",  help='Header for peptide name in peptide metadata file.')
    parser.add_argument('--parentSeqHeader', default="ParentSeqName",  help='Header for parent sequence name in peptide metadata file.')
    parser.add_argument('--enrPepStartStr', default="",  help='Expected string at the beginning of all enriched peptide files.')
    parser.add_argument('--enrPepExtension', default="_enriched.txt",  help='Expected string at the end of all enriched peptide files.')
    parser.add_argument('--pepThresh', default=1, type=int, help='Minimum number of enriched peptides for a p-value to be calculated for a given taxonomic group.')
    parser.add_argument('--basePval', default=0.01, type=float, help='Base p-value threshold to use for reporting.')
    parser.add_argument('--noBonferroni', default=False, action="store_true", help='Use this flag if you do not want to the pvalue threshold to be adjusted using the Bonferroni method.')
    parser.add_argument('--taxMeta', help='Optional metadata file with information about the taxonomic groups to include in the output files.')
    parser.add_argument('--taxMetaKey', default="MultiProteinCluster", help='Header containing taxonomic group name within taxMeta file.')
    parser.add_argument('--taxMetaHeaders', default="Species_List,Species_Prop,Genus_List,Genus_Prop", help='Comma-separated list of headers from the taxMeta file to include in the output.')

    parser.add_argument('rest', nargs=argparse.REMAINDER)
    
    opts = parser.parse_args()
    
    # Read metadata linking protein targets to tax clusters
    prot2taxD = io.fileDictHeader(opts.protMeta, opts.protNameHeader, opts.taxHeader)

    # Read metadata linking peptides to protein targets
    pep2protD = io.fileDictHeader(opts.pepMeta, opts.pepNameHeader, opts.parentSeqHeader)
    
    # If provided, read in taxonomic metadata
    if opts.taxMeta:
        taxDF = pd.read_csv(opts.taxMeta, header=0, index_col=0, sep="\t")
        opts.taxMetaHeaders = opts.taxMetaHeaders.split(",")
    
    # Create output directory and subdirectories, if they doesn't already exist
    io.makeDir(opts.outDir)
    io.makeDir(f"{opts.outDir}/sigTaxa")
    
    
    # For each taxonomic groups, calculate the proportion of the total peptides it represents
    totalPeps = len(pep2protD)
    tax2pepsD = dd(list)
    for pepN, protN in pep2protD.items():
        if protN in prot2taxD:
            tax2pepsD[prot2taxD[protN]].append(pepN)
    
    propLibD = {k:len(v)/totalPeps for k,v in tax2pepsD.items()}
    
    # Grab all of the enriched peptide files and step through them one by one
    for eachF in glob.glob(f"{opts.enrDir}/{opts.enrPepStartStr}*{opts.enrPepExtension}"):
        print(eachF)
        # Read in enriched peptides
        enrPL = io.fileList(eachF, header=False)
        if enrPL != [" "]:
            totalEnrPeps = len(enrPL)
            
            # Generate dictionary with counts of enriched peptides for each taxonomic group
            countD = dd(int)
            
            # Step through each enriched peptide and add a count to the countD
            for pn in enrPL:
                this_prot = pep2protD[pn]
                if this_prot in prot2taxD:
                    countD[prot2taxD[this_prot]]+=1
            
            # Step through each taxonomic groups and calculate probability of observing this number of peptides, given a poisson distribution
            pvD = {}
            for tn, pc in countD.items():
                if pc >= opts.pepThresh:
                    pv = poisson.sf(pc, totalEnrPeps*propLibD[tn])
                    pvD[tn]=pv
            
            # Correct p-value threshold, unless otherwise specified
            if not opts.noBonferroni and len(pvD)>1:
                pval = opts.basePval/len(pvD)
            else:
                pval = opts.basePval
            
            sigNames = [(pv, tn) for tn, pv in pvD.items() if pv<=pval]
            if len(sigNames) > 0:
                with open(f"{opts.outDir}/sigTaxa/{os.path.basename(eachF)[:-4]}.tsv", "w") as fout:
                    if opts.taxMeta:
                        extraHeaders = "\t".join(opts.taxMetaHeaders)
                        fout.write(f"TaxGroup\tNumEnrPeptides\tPoissonProbPvalue\t{extraHeaders}\n")
                    else:
                        fout.write("TaxGroup\tNumEnrPeptides\tPoissonProbPvalue\n")
                    for pv,tn in sorted(sigNames):
                        if opts.taxMeta:
                            extraInfo = "\t".join([str(taxDF[hd][tn]) for hd in opts.taxMetaHeaders])
                            fout.write(f"{tn}\t{countD[tn]}\t{pv}\t{extraInfo}\n")
                        else:
                            fout.write(f"{tn}\t{countD[tn]}\t{pv}\n")
                
    

#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

