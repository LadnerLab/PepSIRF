Learn more or give us feedback
#!/usr/bin/env python

import optparse, os, subprocess, re

#This script reads in a score matrix and outputs sets of "enriched" probes designed from a given taxon

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-m', '--matrix',  help='Score matrix. [None, REQ]')
    p.add_option('-s', '--samps',  help='File with sample names to keep. Can also provide these as arguments [OPT]')
    p.add_option('-o', '--out',  help='Name for output file [None, REQ]')

    opts, args = p.parse_args()
    
    #Read in Score Matrix
    scoreD = parseCounts(opts.matrix)
    
    if opts.samps:
        sampNames = filelist(opts.samps)
    else:
        sampNames = args
    
    subMat = {x:scoreD[x] for x in sampNames}
    writeCounts(subMat, opts.out)

#----------------------End of main()


def filelist(file):
    l=[]
    with open(file, "r") as fin:
        for line in fin:
            l.append(line.strip("\n"))
    return l

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

def writeCounts(cd, outname):
    probeNames = sorted(cd[list(cd.keys())[0]].keys())
    sampNames =  sorted(list(cd.keys()))
    with open(outname, "w") as fout:
        fout.write("Probe\t%s\n" % ("\t".join(sampNames)))
        for each in probeNames:
            fout.write("%s\t%s\n" % (each, "\t".join([str(cd[x][each]) for x in sampNames])))
###------------------------------------->>>>    

if __name__ == "__main__":
    main()
