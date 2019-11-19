#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import numpy as np
import optparse, os, subprocess, re
import itertools as it

#Generate peptide level aligned fasta files based on mapping of probes onto protein-level alignments

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
#    p.add_option('-e', '--enriched',  help='List of peptide names, one per line. [None, REQ]')
    p.add_option('-s', '--scores',  help='Score matrix. Will serve as y-axis in plot. [None, REQ]')
    p.add_option('-n', "--name", help='Sample name, must match a column in the score matrix. [None, REQ]')
    p.add_option('-p', "--probes", help='File containing a list of enriched probes. [None, REQ]')
    p.add_option('-a', '--align',  help='Probe alignment file. [None, REQ]')
    p.add_option('-w', '--win', default=7,  type='int', help='Window size. [7]')
    p.add_option('-t', '--step', default=3,  type='int', help='Step size between windows. [3]')
    p.add_option('--ovlp', default=0.5,  type='float', help='Probe must cover at least this proportion of window to be included. [0.5]')
    p.add_option('-m', '--map',  help='File containing map to link names in enriched list to those in pep fasta. [OPT]')
    p.add_option('--mapOrder', default=1,  type='int', help='Integer indicate the column # of the name labels in the enriched list. Must be 1 or 2 (1-based)[1]')

    opts, args = p.parse_args()

    #Read in name map, if provided
    if opts.map:
        pepMap = readmap(opts.map, order=opts.mapOrder-1)
        revPepMap = {v:k for k,v in pepMap.items()}
        print("Read Name Map")
        
    #Read in score matrix, if provided
    if opts.plot:
        scoreDict = parseCounts(opts.scores)
        print("Read Scores")

    #Read in enriched probes
    pros = filelist(opts.probes)
    print("Read Enriched Probes")

    #Read in probe map
    proAl = readAlignments(opts.align)
    #Make new version with only the enriched probes
    proAl = downSample(proAl, pros, pepMap, revPepMap)
    print("Read Probe Alignment")
    
    
    #Get min and max positions of aligned probes
    minPos, maxPos = getMinMax(proAl)
    
    # Get all species-level taxIDs in the probe alignment
    ids = allIDs(list(proAl.keys()), pepMap, revPepMap)
    
    
    allYs = {x:[] for x in ids}          #To collect maxZ scores, which will be plotted on the y-axis
    xs = []                              #To collect window positions for the x-axis
    
    #Step through each window in the alignment
    for i in range(minPos, maxPos+1-opts.win+1, opts.step):
        xs.append(i)
        thisWin = {x:0 for x in ids}
        inWin = set(range(i,i+opts.win))     #Positions covered by this window
        for p, inf in proAl.items():         #Step through all of the enriched probes
            # If this probe overlaps with the current window
            if len(inf[0].intersection(inWin))/opts.win >= opts.ovlp:
                #If working with the uncoded names
                if "OXX" in p:
                    taxInf = parse_tax(p)
                elif p in pepMap:
                    taxInf = parse_tax(pepMap[p])
                elif p in revPepMap:
                    taxInf = parse_tax(revPepMap[p])
                else:
                    print("Can't find uncoded name for %s" % p)
                if taxInf:
                    if thisWin[taxInf[0]] < scoreDict[opts.name][p]:
                        thisWin[taxInf[0]] = scoreDict[opts.name][p]
                else:
                    print("Can't find tax info for %s" % p)
                    
        for ti, sc in thisWin.items():
            allYs[ti].append(sc)
    
    #Generate and save plot
    plotData(xs, allYs, opts)
    
#----------------------End of main()

def allIDs(pros, pepMap, revPepMap):
    newL=[]
    for p in pros:
        if "OXX" in p:
            taxInf = parse_tax(p)
        elif p in pepMap:
            taxInf = parse_tax(pepMap[p])
        elif p in revPepMap:
            taxInf = parse_tax(revPepMap[p])
        else:
            print("Can't find uncoded name for %s" % p)
        if taxInf:
            newL.append(taxInf[0])

    return list(set(newL))



def downSample(proAl, pros, pepMap, revPepMap):
    newD={}
    for p in pros:
        if p in proAl:
            newD[p] = proAl[p]
        elif p in pepMap:
            newD[p] = proAl[pepMap[p]]
        elif p in revPepMap:
            newD[p] = proAl[revPepMap[p]]
        else:
            print("Can't locate alignment info for %s" % p)
    return newD
    
def getMinMax(proAl):
    mn=9999999
    mx=-9999999
    for p,inf in proAl.items():
        if min(inf[1]) < mn:
            mn=min(inf[1])
        if max(inf[1]) > mx:
            mx=max(inf[1])
    return mn, mx

def plotData(xs, allYs, opts):
    fig,ax = plt.subplots(1,1,figsize=(15,3),facecolor='w')
    
    #Set color palette
    cols = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    
    #One line per species
    counter=0
    for thisID, ys in allYs.items():
        ax.plot(xs, ys, c=cols[counter], label=thisID)

    #Generate figure legend
    ax.legend()

#    ax.set_yticks(range(len(seqs)))
#    ax.set_yticklabels([x.replace("-", " ") for x in seqs])
#    for tick in ax.get_yticklabels():
#        tick.set_fontname("Courier New")
#        tick.set_fontsize(5)

    # Sample labels
#    ax.set_xticks(range(len(samps)))
#    ax.set_xticklabels([x.split("_")[0] for x in samps], rotation = 90)
#    for tick in ax.get_xticklabels():
#        tick.set_fontname("Courier New")
#        tick.set_fontsize(4)
#    plt.tight_layout()
    fig.savefig("%s.pdf" % opts.out,dpi=200,bbox_inches='tight')

def alignSeqs(sD):
    mn = min([min(x[1]) for x in sD.values()])
    mx = max([max(x[1]) for x in sD.values()])
    fullPos = range(mn,mx+1)
    
    names=[]
    seqs=[]
    toSort=[(v[1][0],k) for k,v in sD.items()]
    for st,n in sorted(toSort):
        thisSeqD = {sD[n][1][i]:sD[n][0][i] for i in range(len(sD[n][0]))}
        thisSeq=[]
        for p in fullPos:
            if p in thisSeqD:
                thisSeq.append(str(thisSeqD[p]))
            else:
                thisSeq.append("-")
        names.append(n)
        seqs.append("".join(thisSeq))
        
    #Check for positions that are all gaps
    allGaps=[i for i in range(len(seqs[0])) if set([x[i] for x in seqs]) == set(["-"])]
    #Make new seqs without the gaps
    if allGaps:
        newSeqs = [noGaps(s,allGaps) for s in seqs]
        return names, newSeqs, "%d-%d" % (mn, mx)
    else:
        return names, seqs, "%d-%d" % (mn, mx)

def noGaps(seq, allGaps):
    newS=[]
    for i,s in enumerate(seq):
        if i not in allGaps:
            newS.append(s)
    return "".join(newS)


def clustProbes(eD, proAl):
    ovlp = {x:[] for x in eD}
    for a, b in it.combinations(eD, 2):
        if proAl[a][0].intersection(proAl[b][0]):
            ovlp[a].append(b)
            ovlp[b].append(a)
#        else:
#            print(proAl[a][1])
#            print(proAl[b][1])
    
    clusts = {}
    for s in ovlp:
        this={}
        this = makeClusts(this, s, ovlp)
        c = tuple(sorted(list(this.keys())))
        if c not in clusts:
            clusts[c] = [s]
        else:
            clusts[c].append(s)
    return clusts

def makeClusts(this,s,ovlp):
    if s not in this:
        this[s]=""
        for each in ovlp[s]:
            if each not in this:
                this = makeClusts(this, each, ovlp)
    return this

def readAlignments(align):
    dd={}
    with open(align, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                dd[cols[0]] = (set(range(int(cols[1]), int(cols[2])+1)), [int(x) for x in cols[3].split("~")])
    return dd

def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

def kmers(seq,k):
    out=[]
    for i in range(len(seq)-k+1):
        out.append(seq[i:i+k])
    return set(out)

def filelist(file):
    l=[]
    with open(file, "r") as fin:
        for line in fin:
            l.append(line.strip("\n"))
    return l

def readmap(file, order=0):
    m={}
    with open(file, "r") as fin:
        for line in fin:
            cols=line.strip("\n").split("\t")
            if not order:
                m[cols[0]] = cols[1]
            else:
                m[cols[1]] = cols[0]
    return m



#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def parse_tax(name):
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_id = oxpat.search(name)
    if tax_id:
        species,genus,family = (tax_id.group(2),tax_id.group(3),tax_id.group(4))
        return species,genus,family
    else:
        #print(name)
        return None

def speciesNames(file):
    id2name={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols=line.strip("\n").split("\t")
                id2name[cols[4]]=cols[0]
    return(id2name)

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
