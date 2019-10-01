#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import optparse, re
import numpy as np

#Used to generate plots showing relative enrichment for a set of probes that are most informative for a given taxID

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-d', '--data', help='Count data file. [None, REQ]')
    p.add_option('-s', '--speciesIDs', help='Comma-separated species IDs of interest. [None, REQ]')
    p.add_option('-l', '--linkage', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/speciesDeconvolution/fullTargetsSpeciesK7_wcounts_2019-09-12_coded.tsv",  help='Fasta file with probe/tag sequences. [/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/speciesDeconvolution/fullTargetsSpeciesK7_wcounts_2019-09-12_coded.tsv]')
    p.add_option('-k', '--kDiff', type='int', default=24, help='Minimum additional # of kmers present for focal species in probe, for probe to be included [24]')
    p.add_option('-m', '--min', type='int', default=50, help='Minimum data value to include in color scale [50]')
    p.add_option('-o', '--out', default="informProbes",  help='String to add to output plots. [informProbes]')
    p.add_option('-a', '--align', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/alignments/AlignmentInfo.txt",  help='File with info on probe level alignments. [/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/alignments/AlignmentInfo.txt]')
    p.add_option('--map', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/library_design/oligos/2019-02-27/encoding/unrepSW_9_30_70_design_combined_wControls_2019-09-12.csv_map",  help='Probe name map. [/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/library_design/oligos/2019-02-27/encoding/unrepSW_9_30_70_design_combined_wControls_2019-09-12.csv_map]')
    p.add_option('--log', action='store_true', default=False, help='Use this to plot log scale. [False]')
    p.add_option('--withLabels', action='store_true', default=False, help='Use this flag to include sample labels on y-axis. [False]')
    opts, args = p.parse_args()
    

    #Read in alignment info
    alInfo = {}
    with open(opts.align, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                alInfo[cols[0]] = ["%s/%s" % (cols[1], x) for x in cols[2:]]

    #Generate dict to translate names    
    mapDict=readmap(opts.map, order=1)


    #Read in linkage map
    linkD = readLinkage(opts.linkage)

    #Read in data
    data = parseCounts(opts.data)
    
    #Set color palatte
    palatte = plt.cm.viridis
    palatte.set_under(color='white')

    
    #Step through each taxID of interest
    for thisID in opts.speciesIDs.split(","):
        
        probes = getInfPro(thisID, linkD, opts)               #Pull out informative probes
        #Sort probes based on alignment position, if alignmnet is available
        if thisID in alInfo:
            thisInfo = parseAlInfo(alInfo[thisID])
            toSort = []
            for pr in probes:
                prLong=mapDict[pr]
                for i, inform in  enumerate(list(thisInfo.values())):
                    if prLong in inform:
                        toSort.append((i,inform[prLong][0],pr))
                probes = [x[2] for x in sorted(toSort)]
        
        samps = list(data.keys())                             #Grab sample names
        if opts.log:
            toPlot = [[np.log10(data[x][y]) for y in probes] for x in samps]#Grab data to plot
        else:
            toPlot = [[data[x][y] for y in probes] for x in samps]#Grab data to plot
        
        
        #Generate and save plot
        fig, ax = plt.subplots(figsize=(15,4))
        #Get info for plot

        if opts.log:
            pcm = ax.pcolormesh(toPlot, cmap=palatte, vmin=np.log10(opts.min))  #Plot pcolormesh
        else:
            pcm = ax.pcolormesh(toPlot, cmap=palatte, vmin=opts.min)  #Plot pcolormesh
        #pcm.set_edgecolor('face')

        yticks = [0.5+x for x in list(range(len(samps)))]
        ax.set_yticks(yticks)
        ax.set_ylim(0,len(samps))
        if opts.withLabels:
            ax.set_yticklabels(samps)
        else:
            ax.set_yticklabels([""]*len(yticks))
#        for item in (ax.get_xticklabels()):
#            item.set_fontsize(18)
        for item in (ax.get_yticklabels()):
            item.set_fontsize(8)
        cbar = fig.colorbar(pcm, extend='max')

        if opts.log:
            cbar.ax.set_ylabel('Log Normalized Read Count', fontsize=16)
        else:
            cbar.ax.set_ylabel('Normalized Read Count', fontsize=18)

        ax.set_xlabel("TaxID=%s Required kDiff=%d" % (thisID, opts.kDiff), fontsize=24)
        ax.set_ylabel("Samples", fontsize=24)

        plt.savefig('%s_%d_%s.pdf' % (thisID, opts.kDiff, opts.out),bbox_inches='tight')
        plt.savefig('%s_%d_%s.png' % (thisID, opts.kDiff, opts.out),dpi=300,bbox_inches='tight')


#----------------------End of main()

def getInfPro(thisID, linkD, opts):
    probes=[]
    for p, info in linkD.items():
        if thisID in info:
            if len(info)==1:
                probes.append(p)
            elif info[thisID] == max(info.values()):
                if info[thisID] - max([info[x] for x in info if x != thisID]) >= opts.kDiff:
                     probes.append(p)
    return probes

def readLinkage(linkFile):
    linkD={}
    with open(linkFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                linkD[cols[0]]={}
                for each in cols[1].split(","):
                    if each: 
                        i, c = each.split(":")
                        linkD[cols[0]][i]=int(c)
    return linkD

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

def parseAlInfo(files):
    alInfo={}
    for f in files:
        alPos = parseAlFile(f)
        geneName = f.split("/")[-1].split("_")[2]
        alInfo[geneName] = alPos
    return alInfo

def parseAlFile(file):
    infD={}
    thisMin=1000000
    thisMax=0
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                infD[cols[0]] = [int(cols[1]), int(cols[2]), [int(x) for x in cols[3].split("~")]]
                    
    return infD
























def writeData(samps, z, pos, spID, sp, gene, minDepth, outstr):
    with open("%s_%s_min%d_%s_%s.txt" % (spID, sp, minDepth, gene, outstr), "w") as fout:
        fout.write("Sample\t%s\n" % ("\t".join([str(int(x)) for x in pos])))
        for i,s in enumerate(samps):
            fout.write("%s\t%s\n" % (s, "\t".join([str(x) for x in z[i]])))




###------------------------------------->>>>    

if __name__ == "__main__":
    main()

