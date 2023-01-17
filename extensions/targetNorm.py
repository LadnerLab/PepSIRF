#!/usr/bin/env python

import optparse
import inout as io

# This script can be used to column sum normalize samples using a subset of the peptides in a given assay

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--input',  help='Matrix of raw read counts. [None, REQ]')
    p.add_option('-t', '--targets',  help='File containing the names of the peptides that should be used for the normalization, one peptide per line. [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file, which will contain normalized read count values [None, REQ]')
    opts, args = p.parse_args()
    
    #Read in target peptides
    tarPeps = io.fileList(opts.targets, header=None)
    
    #Read in raw read counts
    inD = parseCounts(opts.input)
    
    #Generate normalized read counts
    normD = {}
    
    for s, dataD in inD.items():
        normD[s] = {}
        tarSum = sum([dataD[p] for p in tarPeps])
        normMultiplier = 1000000/len(dataD)*len(tarPeps)
        for p, v in dataD.items():
            normD[s][p] = v/tarSum*normMultiplier
    
    writeCounts(normD, opts.out)
    
#----------------------End of main()

def writeCounts(cd, outname):
    probeNames = sorted(cd[list(cd.keys())[0]].keys())
    sampNames =  sorted(list(cd.keys()))
    with open(outname, "w") as fout:
        fout.write("Sequence name\t%s\n" % ("\t".join(sampNames)))
        for each in probeNames:
            fout.write("%s\t%s\n" % (each, "\t".join([str(cd[x][each]) for x in sampNames])))

def parseCounts(countFile, delim="\t"):
    counts={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc == 1:
                names = cols[1:]
                for n in names:
                    counts[n]={}
            else:
                for i, count in enumerate(cols[1:]):
                    counts[names[i]][cols[0]] = float(count)
    return counts


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

