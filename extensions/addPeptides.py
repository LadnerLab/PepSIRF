#!/usr/bin/env python

import optparse
from collections import defaultdict

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-p', '--peptides',  help='File containing a list of every peptide that should be in the output, on per row. [None, REQ]')
    p.add_option('-i', '--input',  help='Input matrix. [None, REQ]')
    p.add_option('-o', "--out", help='Name for output file. [None, REQ]')
    p.add_option("--dontReorder", default=False, action="store_true", help='Use this option if you do not care about the order of the rows in the output. This option will require less memory. [OPT]')

    opts, args = p.parse_args()

    #Read in peptide names
    pepD = fileEmptyDict(opts.peptides, col=0, header=False)
    
    with open(opts.out, "w") as fout:
        #If peptides should be reordered, then read in full input matrix
        if opts.dontReorder:
            with open(opts.input, "r") as fin:
                lc=0
                for line in fin:
                    lc+=1
                    cols = line.rstrip("\n").split("\t")
                    header = cols
                    if lc==1:
                        fout.write(line)
                    else:
                        if cols[0] in pepD:
                            fout.write(line)
                            del(pepD[cols[0]])
            for p in pepD:
                fout.write("%s\t%s\n" % (p, "\t".join(["0"]*(len(header)-1))))
        
        else:

            matD, headerL = readMatrixPepKeys(opts.input)
            fout.write("%s\n" % ("\t".join(headerL)))
            for p in pepD:
                if p in matD:
                    fout.write("%s\t%s\n" % (p, "\t".join(matD[p])))
                else:
                    fout.write("%s\t%s\n" % (p, "\t".join(["0"]*(len(headerL)-1))))

#----------------------End of main()

def readMatrixPepKeys(countFile, delim="\t"):
    cmat={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc == 1:
                header = cols
            else:
                cmat[cols[0]] = cols[1:]
    return cmat, header


def fileEmptyDict(file, col=0, delim="\t", header=True):
    l={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                l[thisCols[col]]=""
    return l

###------------------------------------->>>>    

if __name__ == "__main__":
    main()
