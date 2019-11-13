#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import numpy as np
import optparse, os, subprocess, re
import itertools as it

#Generate peptide level aligned fasta files based on mapping of probes onto protein-level alignments

def main():
    usage = '%prog [options] probes1.txt [probes2.txt ...]'
    p = optparse.OptionParser()
#    p.add_option('-e', '--enriched',  help='List of peptide names, one per line. [None, REQ]')
    p.add_option('-p', '--peps',  help='Fasta file peptide sequences. [None, REQ]')
    p.add_option('-a', '--align',  help='Probe alignment file. [None, REQ]')
    p.add_option('-o', "--outDir", default="pepAligned", help='Name for dirctory in which output files will be generated. This Dirctory should NOT exist already. [pepAligned]')
    p.add_option("--plot", help='If this option is used, along with some type of matrix, then heat map plots will be generated for each alignment. [OPT]')
    p.add_option('-m', '--map',  help='File containing map to link names in enriched list to those in pep fasta. [OPT]')
    p.add_option('--mapOrder', default=1,  type='int', help='Integer indicate the column # of the name labels in the enriched list. Must be 1 or 2 (1-based)[1]')
    p.add_option('-s', '--species',  help='File containing link between species taxonomy IDs and full names. [OPT]')

    opts, args = p.parse_args()

    #Check for output directory and create IF it doesn't already exist
    if not os.path.isdir(opts.outDir):
        os.mkdir(opts.outDir)
    
    #Read in peptides
    pepDict = read_fasta_dict_upper(opts.peps)

    #Read in probe map
    proAl = readAlignments(opts.align)

    #Read in name map, if provided
    if opts.map:
        pepMap = readmap(opts.map, order=opts.mapOrder-1)
        revPepMap = {v:k for k,v in pepMap.items()}

    #Read in species ID map, if provided
    if opts.species:
        id2name = speciesNames(opts.species)

    #Read in score matrix, if provided
    if opts.plot:
        scoreDict = parseCounts(opts.plot)
        print("Read Scores")


    #Step through each list of probes provided
    for each in args:
        print(each)
        #Read in peptides of interest
        peps = filelist(each)
        print(len(peps))
        #Generate clusters
        if opts.map:
            clu = clustProbes([pepMap[k] for k in peps if pepMap[k] in proAl], proAl)
        else:
    #        print([k for k in peps if k in proAl])
            clu = clustProbes([k for k in peps if k in proAl], proAl)
        print(len(clu))
        #Check clusters
        for a,b in clu.items():
            if set(a) != set(b):
                print ("Not all starting places resulted in the same cluster!!!", a, sorted(b))
    
        #String for output files
        outStr="%s_%s" % (each.split("/")[-1].split(".")[0], "_".join(opts.align.split("/")[-1].split("_")[:3]))
    
        #Write aligned fasta files
        for c in clu:
            if opts.map:
                names,seqs,bounds = alignSeqs({x:[pepDict[x], proAl[x][1]] for x in c})
            else:
                names,seqs,bounds = alignSeqs({x:[pepDict[x], proAl[x][1]] for x in c})
            write_fasta(names, seqs, "%s/%s_%s.fasta" % (opts.outDir, outStr, bounds))
            
            #Generate heat map
            if opts.plot:
                heatMap([revPepMap[x] for x in names], seqs, {k:{p:s for p,s in v.items() if pepMap[p] in c} for k,v in scoreDict.items()}, "%s/%s_%s.pdf" % (opts.outDir, outStr, bounds), opts)
    
#----------------------End of main()

def heatMap(names, seqs, dta, outname, opts):
    samps = sorted(dta.keys())
    vals = [[dta[x][y] for x in samps] for y in names]
#    print(vals)
    fig,ax = plt.subplots(1,1,figsize=(8,2),facecolor='w')
    
    #Set color palette
    palatte = plt.cm.viridis
    palatte.set_under(color='white')
    
    #Make plot
    ims = ax.imshow(vals, vmin=8)
    cbar = fig.colorbar(ims)
#    ax.set_xlim(-1,len(samps)+1)
#    ax.set_ylim(-1,len(names)+1)

    #Label rows with sequences
    ax.set_yticks(range(len(seqs)))
    ax.set_yticklabels([x.replace("-", " ") for x in seqs])
    for tick in ax.get_yticklabels():
        tick.set_fontname("Courier New")
        tick.set_fontsize(5)

    # Sample labels
    ax.set_xticks(range(len(samps)))
    ax.set_xticklabels([x.split("_")[0] for x in samps], rotation = 90)
    for tick in ax.get_xticklabels():
#        tick.set_fontname("Courier New")
        tick.set_fontsize(4)
    plt.tight_layout()
    fig.savefig(outname,dpi=200,bbox_inches='tight')

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



def mergeClusts(clusts, k, indivThresh, clustThresh):
    pairs2merge = []
    for i in range(len(clusts)):
        for j in range(i+1, len(clusts)):
            hits=0
            comps=0
            for n1,s1 in clusts[i].items():
                for n2,s2 in clusts[j].items():
                    comps+=1
                    if kmerOvlp(s1, s2, k)>=indivThresh:
                        hits+=1
            if hits/comps >= clustThresh:
                pairs2merge.append([i,j])
    
    groups2merge = combPairs(pairs2merge)
    
    newClusts = []
    singles = []
    merged = []
    for each in groups2merge:
        newClusts.append({})
        for a in each:
            merged.append(a)
            for k,v in clusts[a].items():
                newClusts[-1][k] = v
    for i,info in enumerate(clusts):
        if i not in merged:
            if len(info) == 1:
                singles.append(info)
            else:
                newClusts.append(info)
    return newClusts, singles

def combPairs(pairs):
    d={}
    groups = pairs[::]
    for each in pairs:
        for x in each:
            d[x] = d.get(x, 0) + 1
    for each, count in d.items():
        if count>1:
            groups = combine(each, groups)
    return groups
            
def combine(focal, gList):
    newL = []
    combo = []
    for each in gList:
        if focal in each:
            combo+=each
        else:
            newL.append(each)
    newL.append(list(set(combo)))
    return newL


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
        return family,genus,species
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

