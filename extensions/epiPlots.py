#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import numpy as np
import optparse, os, subprocess, re
import itertools as it
import inout as io
from collections import defaultdict
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#Generate epitope-level heat maps for some provided score

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
#    p.add_option('-e', '--enriched',  help='List of peptide names, one per line. [None, REQ]')
    p.add_option('-i', '--inputInfo',  help='Tab-delimited file indicating the sample names and enriched probe lists. If you want replicates to be averaged, probe name as a unique starting string common to all replicates. [None, REQ]')
    p.add_option('-m', '--meta',  help='Metadata file. [None, REQ]')
    p.add_option('-o', "--outDir", default="pepAligned", help='Name for dirctory in which output files will be generated. This Dirctory should NOT exist already. [pepAligned]')
    p.add_option("--plot", help='If this option is used, along with some type of matrix, then heat map plots will be generated for each alignment. [None, REQ]')
    p.add_option("--minToPlot", default=1, type="float", help='Minimum score to include in plots. [1]')
    p.add_option("--plotLog", default=False, action="store_true", help='Use this flag if you want the y-axis for plots to be on the log scale. [False]')
    p.add_option("--addNegNorm", help='Use this flag if you want to add info about the norm read counts for negative controls y-axis labels. Argument should be tab-delimited file. First column should be comma-sep list of negative control names, the second the name of the norm read count matrix. [None, OPT]')

    opts, args = p.parse_args()

    #Check for output directory and create IF it doesn't already exist
    if not os.path.isdir(opts.outDir):
        os.mkdir(opts.outDir)
    
    #Get neg norm counts, if provided:
    normNeg = {}
    if opts.addNegNorm:
        with open(opts.addNegNorm, "r") as fin:
            lc=0
            for line in fin:
                lc+=1
                if lc==1:
                    cols = line.rstrip("\n").split("\t")
                    negNames = cols[0].split(",")
                    normDict = parseCounts(cols[1])
        for p in normDict[negNames[0]]:
            normNeg[p] = np.mean([normDict[x][p] for x in negNames])
    
    #Read in info from metadata file
    pepD = io.fileDictHeader(opts.meta, "CodeName", "Peptide")
    catD=io.fileDictHeader(opts.meta, "CodeName", "Category")
    protD=io.fileDictHeader(opts.meta, "CodeName", "Protein")
    startD=io.fileDictHeader(opts.meta, "CodeName", "AlignStart")
    startD={k:int(v) for k,v in startD.items() if v}
    stopD=io.fileDictHeader(opts.meta, "CodeName", "AlignStop")
    stopD={k:int(v) for k,v in stopD.items() if v}

    #Read in score matrix for generating heat maps, if provided
    scoreDict = parseCounts(opts.plot)

    #Reads in names and probe lists
    enrichD = io.fileDictHeader(opts.inputInfo, "Sample", "PepFile")

    #Generate replicate lists for each name provided
    repD = {s:[n for n in scoreDict if n.startswith(s)] for s in enrichD}
    
    #Step through each list of probes provided
    masterByProt = defaultdict(list)    #To hold info
    for sName, each in enrichD.items():
        #Read in peptides of interest
        peps = filelist(each)
        
        for p in peps:
            if catD[p] != "Control":
                masterByProt[protD[p]].append(p)
        
    for prot, pL in masterByProt.items():
        clu = clustProbes(pL, pepD, startD, stopD)
        #Check clusters
        for a,b in clu.items():
            if set(a) != set(b):
                print ("Not all starting places resulted in the same cluster!!!", a, sorted(b))

        #String for output files
        outStr="%s_master" % (prot)

        for c in clu:
            names,seqs,bounds = alignSeqs({x:[pepD[x], list(range(startD[x], stopD[x]))] for x in c})

            #Write aligned fasta files
#            write_fasta(names, seqs, "%s/%s_%s.fasta" % (opts.master, outStr, bounds))

            #Generate heat map
            dta = {sN:{p:0.0001 for p in c} for sN in repD}
            for sN, each in enrichD.items():
                #Read in peptides of interest
                peps = filelist(each)
        
                for p in peps:
                    if catD[p] != "Control":
                        avgV = np.mean([scoreDict[x][p] for x in repD[sN]])
                        if avgV == float('inf'):
                            avgV = 10000
                        dta[sN][p] = avgV
            
            if opts.addNegNorm:
                seqsPlus = ["%.0f %s" % (normNeg[names[i]], seqs[i]) for i in range(len(seqs))]
                heatMap(names, seqsPlus, dta, "%s/%s_%s.pdf" % (opts.outDir, outStr, bounds), opts.minToPlot, opts.plotLog)
            else:
                heatMap(names, seqs, dta, "%s/%s_%s.pdf" % (opts.outDir, outStr, bounds), opts.minToPlot, opts.plotLog)


    
#----------------------End of main()

def heatMap(names, seqs, dta, outname, minToPlot, useLog):
    samps = sorted(dta.keys())
    vals = [[dta[x][y] for x in samps] for y in names]
    if useLog:
        vals = [[np.log10(y) for y in x] for x in vals]
        minToPlot = np.log10(minToPlot)
#    print(vals)
    w = len(seqs[0])*0.25+len(samps)*0.2
    h = len(seqs)*0.5
    print(outname, w,h)
    fig,ax = plt.subplots(1,1,figsize=(w,h),facecolor='w')
    
    #Set color palette
    palatte = plt.cm.viridis
    palatte.set_under(color='white')
    
    #Make plot
    ims = ax.imshow(vals, vmin=minToPlot)
    
    
    # Make legend
    axins = inset_axes(ax,
                   width="10%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
    cbar = fig.colorbar(ims, cax=axins)
#    ax.set_xlim(-1,len(samps)+1)
#    ax.set_ylim(-1,len(names)+1)

    #Label rows with sequences
    ax.set_yticks(range(len(seqs)))
    ax.set_yticklabels([x.replace("-", " ") for x in seqs])
    for tick in ax.get_yticklabels():
        tick.set_fontname("Courier New")
        tick.set_fontsize(14)

    # Sample labels
    ax.set_xticks(range(len(samps)))
    ax.set_xticklabels([x.split("_")[0] for x in samps], rotation = 90)
    for tick in ax.get_xticklabels():
#        tick.set_fontname("Courier New")
        tick.set_fontsize(14)
#    plt.tight_layout()
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

def clustProbes(focal, pepD, startD, stopD):
#def clustProbes(eD, proAl):
    ovlp = {x:[] for x in focal}
    for a, b in it.combinations(focal, 2):
        aPos = set(range(startD[a], stopD[a]))
        bPos = set(range(startD[b], stopD[b]))
        if aPos.intersection(bPos):
            ovlp[a].append(b)
            ovlp[b].append(a)
    
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
