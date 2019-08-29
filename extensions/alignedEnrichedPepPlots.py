#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

import optparse, re
import numpy as np

#Used to generate plots showing the spatial orientation of enriched probes across an an alignment for multiple individuals and viral species

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-d', '--data', help='Count data file. [None, REQ]')
    p.add_option('-s', '--speciesIDs', help='Comma-separated species IDs of interest. [None, REQ]')
    p.add_option('-p', '--probes', default="/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/PV1_10K3000_53_encoded.fna",  help='Fasta file with probe/tag sequences. [/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/PV1_10K3000_53_encoded.fna]')
    p.add_option('-m', '--map', default="/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/unrepSW_9_30_70_design_combined_wControls.csv_map",  help='Probe name map. [/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/unrepSW_9_30_70_design_combined_wControls.csv_map]')
    p.add_option('-t', '--taxa', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Manuscripts/PV-PepSeq_Design-Testing/Tables/speciesIDs_2019-06-21.txt", help='Taxa info to link IDs to names [/Volumes/GoogleDrive/Shared drives/LadnerLab/Manuscripts/PV-PepSeq_Design-Testing/Tables/speciesIDs_2019-06-21.txt]')
    p.add_option('-a', '--alignInfo', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/alignments/AlignmentInfo.txt", help='Contains info for species-level seq alignments. [/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/alignments/AlignmentInfo.txt]')
    p.add_option('-n', '--minDepth', type='float', default=100, help='Minimum count to be included in plot. [100]')
    p.add_option('--readDepthToo', action='store_true', default=False, help='Use this flag to also generate a plot showing combined read depth at each position. This is most appropriate with normalized read count data. [False]')
    opts, args = p.parse_args()
    

    #Read in probes
    names,tseqs = read_fasta_lists(opts.probes)
    #Creat dict with seqs as keys, names as values
    tagnames={tseqs[i]:names[i] for i in range(len(names))}
    
    #Generate dict to translate names    
    mapDict=readmap(opts.map, order=1)

    #Read in info that will allow taxID to be linked to a species name
    id2name={}
    with open(opts.taxa, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols=line.strip("\n").split("\t")
                id2name[cols[4]]=cols[0]

    #Read in data
    data = parseCounts(opts.data)
    
    #Read in alignment data
    alInfo = {}
    with open(opts.alignInfo, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                alInfo[cols[0]] = ["%s/%s" % (cols[1], x) for x in cols[2:]]

    
    #Make probe count plots
    #By default, any sample that starts with "Super" is excluded. These are expected to be negative controls
    for each in opts.speciesIDs.split(","):
        plotAlignHits([x for x in list(data.keys()) if not x.startswith("Super")], each, alInfo, data, mapDict, id2name, opts.minDepth)
    

    #Make read count depth plots (will be most appropriate with normalized read counts)
    #By default, any sample that starts with "Super" is excluded. These are expected to be negative controls
    if opts.readDepthToo:
        for each in opts.speciesIDs.split(","):
            plotAlignDepth([x for x in list(data.keys()) if not x.startswith("Super")], each, alInfo, data, mapDict, id2name, opts.minDepth)

#----------------------End of main()

def plotAlignDepth(samps, spID, alInfo, data, mapDict, id2name, minDepth, alFasta=False):
    proNames = speciesProbes(mapDict, spID)
    thisInfo = parseAlInfo(alInfo[spID])

    #Set color palatte
    palatte = plt.cm.viridis
    palatte.set_under(color='white')

    #If there is only one alignment
    if len(thisInfo) == 1:
        fig, ax = plt.subplots(figsize=(20,4))
        #Get info for plot
        x,y,z = plotInfoDepth(thisInfo[0], samps, mapDict, data, proNames, minDepth)
        
        #Plot pcolormesh
        pcm = ax.pcolormesh(x, y, z, cmap=palatte, vmin=0.01)
        #pcm.set_edgecolor('face')

        yticks = [0.5+x for x in list(range(1, len(samps)+1))]
        ax.set_yticks(yticks)
        ax.set_ylim(1,len(samps)+1)
        ax.set_yticklabels([""]*len(yticks))
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(18)
        cbar = fig.colorbar(pcm, extend='max')
        cbar.ax.set_ylabel('Normalized Read Count\n(min=%d)' % minDepth, fontsize=18)

        #ax.set_xlabel("Generation Time\n(days between infection of index and secondary cases)", fontsize=18)
        ax.set_xlabel("%s%s" % (thisInfo[0][0][0].upper(), thisInfo[0][0][1:]), fontsize=24)
        ax.set_ylabel("Individuals", fontsize=24)
        
        #Write out data used to make plot
        writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[0][0], minDepth, "depthData")

    else:
        geneSizes = [x[3]-x[2]+1 for x in thisInfo]
        fig, ax = plt.subplots(1, len(thisInfo), figsize=(20,4), sharey=True, gridspec_kw={'width_ratios': geneSizes, 'wspace':0.05})
        
        plotData=[]
        for each in thisInfo:
            #Get info for plot
            plotData.append(plotInfoDepth(each, samps, mapDict, data, proNames, minDepth))
        
        maxZ=0
        for each in plotData:
            for every in each[2]:
                if max(every)>maxZ:
                    maxZ = max(every)
        
        for i, each in enumerate(plotData):
            x,y,z = each
            #Plot pcolormesh
            pcm = ax[i].pcolormesh(x, y, z, cmap=palatte, vmin=0.01, vmax=maxZ)

            yticks = [0.5+x for x in list(range(1, len(samps)+1))]
            ax[i].set_yticks(yticks)
            ax[i].set_ylim(1,len(samps)+1)
            ax[i].set_yticklabels([""]*len(yticks))
            xticks = [0] + [a for a in x[0] if a%500==0]
            if len(xticks)==1:
                xticks.append(max(x[0]))
            ax[i].set_xticks(xticks)


#            for item in (ax[i].get_xticklabels()):
#                item.set_fontsize(18)
            ax[i].set_xlabel("%s%s" % (thisInfo[i][0][0].upper(), thisInfo[i][0][1:]), fontsize=16)

            if i==0:
                ax[i].set_ylabel("Individuals", fontsize=24)
            elif i==len(thisInfo)-1:
                cbar = fig.colorbar(pcm, extend='max', ax=ax)
                cbar.ax.set_ylabel('Normalized Read Count\n(min=%d)' % minDepth, fontsize=18)

            #Write out data used to make plot
            writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[i][0], minDepth, "depthData")
            

    plt.suptitle("%s" % id2name[spID], fontsize=30)
#    plt.tight_layout()
    plt.savefig('%s_%s_min%d_Depth.pdf' % (spID, id2name[spID].replace(" ", "-"), minDepth),dpi=200,bbox_inches='tight')
    plt.savefig('%s_%s_min%d_Depth.png' % (spID, id2name[spID].replace(" ", "-"), minDepth),dpi=200,bbox_inches='tight')


def plotAlignHits(samps, spID, alInfo, data, mapDict, id2name, minDepth, alFasta=False):
    proNames = speciesProbes(mapDict, spID)
    thisInfo = parseAlInfo(alInfo[spID])

    #Set color palatte
    palatte = plt.cm.plasma
    palatte.set_under(color='white')

    #If there is only one alignment
    if len(thisInfo) == 1:
        fig, ax = plt.subplots(figsize=(20,4))
        #Get info for plot
        x,y,z = plotInfo(thisInfo[0], samps, mapDict, data, proNames, minDepth)
        
        #Plot pcolormesh
        pcm = ax.pcolormesh(x, y, z, cmap=palatte, vmin=0.01)
        #pcm.set_edgecolor('face')

        yticks = [0.5+x for x in list(range(1, len(samps)+1))]
        ax.set_yticks(yticks)
        ax.set_ylim(1,len(samps)+1)
        ax.set_yticklabels([""]*len(yticks))
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(18)
        cbar = fig.colorbar(pcm, extend='max')
        cbar.ax.set_ylabel('# Enriched Probes\n(min=%d)' % minDepth, fontsize=18)

        #ax.set_xlabel("Generation Time\n(days between infection of index and secondary cases)", fontsize=18)
        ax.set_xlabel("%s%s" % (thisInfo[0][0][0].upper(), thisInfo[0][0][1:]), fontsize=24)
        ax.set_ylabel("Individuals", fontsize=24)
        
        #Write out data used to make plot
        writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[0][0], minDepth, "data")

    else:
        geneSizes = [x[3]-x[2]+1 for x in thisInfo]
        fig, ax = plt.subplots(1, len(thisInfo), figsize=(20,4), sharey=True, gridspec_kw={'width_ratios': geneSizes, 'wspace':0.05})
        
        plotData=[]
        for each in thisInfo:
            #Get info for plot
            plotData.append(plotInfo(each, samps, mapDict, data, proNames, minDepth))
        
        maxZ=0
        for each in plotData:
            for every in each[2]:
                if max(every)>maxZ:
                    maxZ = max(every)
        
        for i, each in enumerate(plotData):
            x,y,z = each
            #Plot pcolormesh
            pcm = ax[i].pcolormesh(x, y, z, cmap=palatte, vmin=0.01, vmax=maxZ)

            yticks = [0.5+x for x in list(range(1, len(samps)+1))]
            ax[i].set_yticks(yticks)
            ax[i].set_ylim(1,len(samps)+1)
            ax[i].set_yticklabels([""]*len(yticks))
            xticks = [0] + [a for a in x[0] if a%500==0]
            if len(xticks)==1:
                xticks.append(max(x[0]))
            ax[i].set_xticks(xticks)


#            for item in (ax[i].get_xticklabels()):
#                item.set_fontsize(18)
            ax[i].set_xlabel("%s%s" % (thisInfo[i][0][0].upper(), thisInfo[i][0][1:]), fontsize=16)

            if i==0:
                ax[i].set_ylabel("Individuals", fontsize=24)
            elif i==len(thisInfo)-1:
                cbar = fig.colorbar(pcm, extend='max', ax=ax)
                cbar.ax.set_ylabel('# Enriched Probes\n(min=%d)' % minDepth, fontsize=18)

            #Write out data used to make plot
            writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[i][0], minDepth, "data")
            

    plt.suptitle("%s" % id2name[spID], fontsize=30)
#    plt.tight_layout()
    plt.savefig('%s_%s_min%d.pdf' % (spID, id2name[spID].replace(" ", "-"), minDepth),dpi=200,bbox_inches='tight')
    plt.savefig('%s_%s_min%d.png' % (spID, id2name[spID].replace(" ", "-"), minDepth),dpi=200,bbox_inches='tight')

def writeData(samps, z, pos, spID, sp, gene, minDepth, outstr):
    with open("%s_%s_min%d_%s_%s.txt" % (spID, sp, minDepth, gene, outstr), "w") as fout:
        fout.write("Sample\t%s\n" % ("\t".join([str(int(x)) for x in pos])))
        for i,s in enumerate(samps):
            fout.write("%s\t%s\n" % (s, "\t".join([str(x) for x in z[i]])))

def plotInfoDepth(info, samps, mapDict, data, proNames, minDepth):
    geneName, alPos, thisMin, thisMax = info
    covPos = range(thisMin, thisMax+1)
    posCounts = {n:{x:0 for x in covPos} for n in samps}

    #Counts for replicates will be averaged
    for each in samps:
        names=[x for x in data.keys() if x.startswith(each)]
        for p in proNames:
            if mapDict[p] in alPos:
                avgCount = np.mean([data[x][p] for x in names])
                if avgCount >= minDepth:
                    for pos in alPos[mapDict[p]][2]:
                        posCounts[each][pos]+=avgCount

    posDepthLists = []
    for s in samps:
        posDepthLists.append([posCounts[s][x] for x in covPos])
                    
    xvec = np.linspace(min(covPos),max(covPos)+1,num=len(covPos)+1)
    yvec = np.linspace(1,len(samps)+1,num=len(samps)+1)

    xmat = np.array([xvec]*(len(samps)+1))
    ymat = np.column_stack([yvec[::-1]]*(len(covPos)+1))

    z = np.array(posDepthLists)
    return xmat, ymat, z

def plotInfo(info, samps, mapDict, data, proNames, minDepth):
    geneName, alPos, thisMin, thisMax = info
    covPos = range(thisMin, thisMax+1)
    posCounts = {n:{x:0 for x in covPos} for n in samps}

    #Counts for replicates will be averaged
    for each in samps:
        names=[x for x in data.keys() if x.startswith(each)]
        for p in proNames:
            if mapDict[p] in alPos:
                avgCount = np.mean([data[x][p] for x in names])
                if avgCount >= minDepth:
                    for pos in alPos[mapDict[p]][2]:
                        posCounts[each][pos]+=1

    posCountsLists = []
    for s in samps:
        posCountsLists.append([posCounts[s][x] for x in covPos])
                    
    xvec = np.linspace(min(covPos),max(covPos)+1,num=len(covPos)+1)
    yvec = np.linspace(1,len(samps)+1,num=len(samps)+1)

    xmat = np.array([xvec]*(len(samps)+1))
    ymat = np.column_stack([yvec[::-1]]*(len(covPos)+1))

    z = np.array(posCountsLists)
    return xmat, ymat, z



def fileList(file):
    l=[]
    with open(file, "r") as fin:
        for line in fin:
            l.append(line.rstrip("\n"))
    return(l)


def parse_tax(name):
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_id = oxpat.search(name)
    if tax_id:
        species,genus,family = (tax_id.group(2),tax_id.group(3),tax_id.group(4))
        return family,genus,species
    else:
        #print(name)
        return None

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

def speciesProbes(mapDict, spID):
    spProbes=[]
    for short, long in mapDict.items():
        taxInfo = parse_tax(long)
        if taxInfo:
            if taxInfo[2] == spID:
                spProbes.append((int(long.split("_")[-2]), short))
    return [x[1] for x in sorted(spProbes)]

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
    alInfo=[]
    for f in files:
        alPos, thisMin, thisMax = parseAlFile(f)
        geneName = f.split("/")[-1].split("_")[2]
        alInfo.append([geneName, alPos, thisMin, thisMax])
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
                if int(cols[1]) < thisMin: 
                    thisMin=int(cols[1])
                if int(cols[2]) > thisMax: 
                    thisMax=int(cols[2])
                    
    return infD, thisMin, thisMax


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

