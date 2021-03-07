#!/usr/bin/env python

import argparse, glob, os, sys
import inout as io             #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script reads in one or more lists of peptides, along with a metadata file 
# And generates a matrix of counts, one column for each species from which the peptides were designed and one row for each input peptide lists

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pepFile", help="A file containing a list of peptide names", nargs='*')
    parser.add_argument("-n", "--name", help="Column name for peptide names in metadata file", default="CodeName")
    parser.add_argument("-c", "--cat", help="Column name for category of interest in metadata file", default="Species")
    
    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-o", "--out", help="Output matrix file name", required=True)
    reqArgs.add_argument("-m", "--meta", help="Metadata file name", required=True)

    args = parser.parse_args()

    # Read in info from metadata file
    catD = io.fileDictHeader(args.meta, args.name, args.cat)
    
    # Create dictionary to hold counts
    countD = defaultdict(dict)
    
    # Step through each peptides file
    for pF in args.pepFile:
        # Read in peptide names
        peps = io.fileList(pF, header=False)
        for p in peps:
            c = catD[p]
            countD[c][pF] = countD[c].get(pF, 0) + 1
    
    # Sorted list of categories
    allCats = sorted(list(countD.keys()))
    
    #Write output file
    with open(args.out, "w") as fout:
        fout.write("File\t%s\n" % ("\t".join(allCats)))
        for pF in args.pepFile:
            counts = [str(countD[c][pF]) if pF in countD[c] else "0" for c in allCats]
            fout.write("%s\t%s\n" % (pF, "\t".join(counts)))

#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

