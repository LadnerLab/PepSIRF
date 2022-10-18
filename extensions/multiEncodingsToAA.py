#!/usr/bin/env python

import argparse
import numpy as np

def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg_parser.add_argument('--delim', default = "-", help = "The portion of the encoding names before this character should correspond to the peptide-level name")

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-m', '--matrix', help = "Input count matrix within multiple encodings per peptide.", required=True )
    reqArgs.add_argument( '-o', '--out', help = "Name for output matrix with peptide-level counts", required=True )

    args = arg_parser.parse_args()
    
    # Dictionary to hold and sum counts
    outD = {}
    
    with open(args.out, "w") as fout:
    #Step through input file line by line
        lc=0
        with open(args.matrix, "r") as fin:
            for line in fin:
                lc+=1
                # Header line
                if lc == 1:
                    fout.write(line)
                else:
                    cols = line.rstrip("\n").split("\t")
                    pep = cols[0].split(args.delim)[0]
                    if pep not in outD:
                        outD[pep] = np.array([float(x) for x in cols[1:]])
                    else:
                        outD[pep] = outD[pep] + np.array([float(x) for x in cols[1:]])
        #Write out new counts
        for pName, counts in outD.items():
            fout.write("%s\t%s\n" % (pName, "\t".join([str(x) for x in counts])))
        
#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

