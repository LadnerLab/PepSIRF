#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import argparse, re
import inout as io
import numpy as np

#Used to generate plots showing the spatial orientation of enriched probes across an an alignment for multiple individuals and viral species

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pepFile", help="A file containing a list of peptide names", nargs='*')
    parser.add_argument("--name", help="Column name for peptide names in metadata file", default="CodeName")
    parser.add_argument("--sid", help="Column name for species IDs in metadata file", default="SpeciesID")
    parser.add_argument("--species", help="Column name for species names in metadata file", default="Species")
    parser.add_argument('-a', '--alignInfo', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/analysis/alignments/AlignmentInfoCoded.txt", help='Contains info for species-level seq alignments.')
    parser.add_argument('--sampleOrder', help='Can provide a plain text file with the order you would like the samples to appear in the plot.')
    parser.add_argument('--annots', help='File with annotation info. If provided, a graphical representation of the different proteins with also be plotted.')
    
    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-o", "--out", help="Output matrix file name", required=True)
    reqArgs.add_argument("-m", "--meta", help="Metadata file name", required=True)
    reqArgs.add_argument('-s', '--speciesIDs', help='Comma-separated species IDs of interest.', required=True)

    args = parser.parse_args()

#    p.add_option('-p', '--probes', default="/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/PV1_10K3000_53_encoded.fna",  help='Fasta file with probe/tag sequences. [/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/PV1_10K3000_53_encoded.fna]')
#    p.add_option('-m', '--map', default="/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/unrepSW_9_30_70_design_combined_wControls.csv_map",  help='Probe name map. [/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/encoding/unrepSW_9_30_70_design_combined_wControls.csv_map]')
#    p.add_option('-t', '--taxa', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Manuscripts/PV-PepSeq_Design-Testing/Tables/speciesIDs_2019-06-21.txt", help='Taxa info to link IDs to names [/Volumes/GoogleDrive/Shared drives/LadnerLab/Manuscripts/PV-PepSeq_Design-Testing/Tables/speciesIDs_2019-06-21.txt]')
#    p.add_option('--withLabels', action='store_true', default=False, help='Use this flag to include sample labels on y-axis. [False]')


    # Read in info from metadata file
    sidD = io.fileDictHeader(args.meta, args.name, args.sid)
    id2name = io.fileDictHeader(args.meta, args.sid, args.species)

    #Read in sample order for plots, IF provided
    sampOrderList = io.fileList(args.sampleOrder, header=False)

    #Read in probes
#    names,tseqs = read_fasta_lists(opts.probes)
    #Creat dict with seqs as keys, names as values
#    tagnames={tseqs[i]:names[i] for i in range(len(names))}
    
    #Generate dict to translate names    
#    mapDict=readmap(opts.map, order=1)

    #Read in alignment data
    alInfo = {}
    with open(args.alignInfo, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols = line.rstrip("\n").split("\t")
                alInfo[cols[0]] = ["%s/%s" % (cols[1], x) for x in cols[2:]]

    #Read in annotation info, if provided
    annotD = {}
    if args.annots:
        with open(args.annots, "r") as fin:
            lc=0
            for line in fin:
                lc+=1
                if lc>1:
                    cols = line.rstrip("\n").split("\t")
                    if cols[0] not in annotD:
                        annotD[cols[0]] = {}
                    annotD[cols[0]][cols[1]] = [x.split(",") for x in cols[2].split("~")]
    
    #Make probe count plots
    #By default, any sample that starts with "Super" is excluded. These are expected to be negative controls
    for sid in args.speciesIDs.split(","):
        if args.sampleOrder:
            plotAlignHits(sampOrderList, sid, alInfo, id2name[sid], annotD, args)
        else:
            plotAlignHits([x for x in list(data.keys()) if not x.startswith("Super")], each, alInfo, data, mapDict, id2name, opts.minDepth, annotD, opts)
    

#----------------------End of main()

def plotAnnots(aD, ax):
    prevY=0
    prevStop=0
    for name, genes in aD.items():
        for start, end, g in genes:
            if int(start)<prevStop:
                y = prevY+0.5
            else:
                y= 0

            p = patches.Rectangle((int(start), y), int(end)-int(start)+1, 1,fill=False, clip_on=False)
            ax.add_patch(p)
            
            ax.text(0.5*(int(start)+int(end)), 0.5*(y+y+1), g,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=15, color='black')
        
        ax.set_axis_off()



def plotAlignHits(samps, spID, alInfo, sname, annotD, opts, alFasta=False):
    
    #Alignment info for this species
    thisInfo = parseAlInfo(alInfo[spID])

    #Set color palatte
    palatte = plt.cm.plasma
    palatte.set_under(color='white')

    #If there is only one alignment
    if len(thisInfo) == 1:
        #Get info for plot
        x,y,z = plotInfo(thisInfo[0], samps)

        if spID not in annotD:
            fig, ax = plt.subplots(figsize=(20,4))
        
            #Plot pcolormesh
            pcm = ax.pcolormesh(x, y, z, cmap=palatte, vmin=0.01)

            yticks = [0.5+x for x in list(range(1, len(samps)+1))]
            ax.set_yticks(yticks)
            ax.set_ylim(1,len(samps)+1)
#            if opts.withLabels:
#                ax.set_yticklabels(samps[::-1])
#            else:
            ax.set_yticklabels([""]*len(yticks))
            for item in (ax.get_xticklabels()):
                item.set_fontsize(18)
            for item in (ax.get_yticklabels()):
                item.set_fontsize(10)
            cbar = fig.colorbar(pcm, extend='max')
            cbar.ax.set_ylabel('# Enriched Peptides', fontsize=18)

            #ax.set_xlabel("Generation Time\n(days between infection of index and secondary cases)", fontsize=18)
#            ax.set_xlabel("%s%s" % (thisInfo[0][0][0].upper(), thisInfo[0][0][1:]), fontsize=24)
            ax.set_ylabel("Individuals", fontsize=24)
        
        else:
            fig, ax = plt.subplots(2,1, figsize=(20,4), sharex=True, gridspec_kw={'height_ratios': [1,5], 'hspace':0.1})
            plotAnnots(annotD[spID], ax[0])
        
            #Plot pcolormesh
            pcm = ax[1].pcolormesh(x, y, z, cmap=palatte, vmin=0.01)

            yticks = [0.5+x for x in list(range(1, len(samps)+1))]
            ax[1].set_yticks(yticks)
            ax[1].set_ylim(1,len(samps)+1)
#            if opts.withLabels:
#                ax[1].set_yticklabels(samps[::-1])
#            else:
            ax[1].set_yticklabels([""]*len(yticks))
            for item in (ax[1].get_xticklabels()):
                item.set_fontsize(18)
            for item in (ax[1].get_yticklabels()):
                item.set_fontsize(10)
            plt.xlim(0,len(x[0]))
            cbar = fig.colorbar(pcm, extend='max', ax=ax)
            cbar.ax.set_ylabel('# Enriched Peptides', fontsize=18)

            #ax.set_xlabel("Generation Time\n(days between infection of index and secondary cases)", fontsize=18)
#            ax[1].set_xlabel("%s%s" % (thisInfo[0][0][0].upper(), thisInfo[0][0][1:]), fontsize=24)
            ax[1].set_ylabel("Individuals", fontsize=24)

        #Write out data used to make plot
#        writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[0][0], minDepth, "data")
# 
#     else:
#         geneSizes = [x[3]-x[2]+1 for x in thisInfo]
#         fig, ax = plt.subplots(1, len(thisInfo), figsize=(20,4), sharey=True, gridspec_kw={'width_ratios': geneSizes, 'wspace':0.05})
#         
#         plotData=[]
#         for each in thisInfo:
#             #Get info for plot
#             plotData.append(plotInfo(each, samps, mapDict, data, proNames, minDepth))
#         
#         maxZ=0
#         for each in plotData:
#             for every in each[2]:
#                 if max(every)>maxZ:
#                     maxZ = max(every)
#         
#         for i, each in enumerate(plotData):
#             x,y,z = each
#             #Plot pcolormesh
#             pcm = ax[i].pcolormesh(x, y, z, cmap=palatte, vmin=0.01, vmax=maxZ)
# 
#             yticks = [0.5+x for x in list(range(1, len(samps)+1))]
#             ax[i].set_yticks(yticks)
#             ax[i].set_ylim(1,len(samps)+1)
#             ax[i].set_yticklabels([""]*len(yticks))
#             xticks = [0] + [a for a in x[0] if a%500==0]
#             if len(xticks)==1:
#                 xticks.append(max(x[0]))
#             ax[i].set_xticks(xticks)
# 
# 
# #            for item in (ax[i].get_xticklabels()):
# #                item.set_fontsize(18)
#             ax[i].set_xlabel("%s%s" % (thisInfo[i][0][0].upper(), thisInfo[i][0][1:]), fontsize=16)
# 
#             if i==0:
#                 ax[i].set_ylabel("Individuals", fontsize=24)
#             elif i==len(thisInfo)-1:
#                 cbar = fig.colorbar(pcm, extend='max', ax=ax)
#                 cbar.ax.set_ylabel('# Enriched Probes\n(min=%d)' % minDepth, fontsize=18)
# 
#             #Write out data used to make plot
#             writeData(samps, z, x[0][:-1], spID, id2name[spID], thisInfo[i][0], minDepth, "data")
            

    plt.suptitle("%s" % sname, fontsize=30)
#    plt.tight_layout()
    plt.savefig('%s_%s.pdf' % (spID, sname.replace(" ", "-")),dpi=200,bbox_inches='tight')
    plt.savefig('%s_%s.png' % (spID, sname.replace(" ", "-")),dpi=200,bbox_inches='tight')

def writeData(samps, z, pos, spID, sp, gene, minDepth, outstr):
    with open("%s_%s_min%d_%s_%s.txt" % (spID, sp, minDepth, gene, outstr), "w") as fout:
        fout.write("Sample\t%s\n" % ("\t".join([str(int(x)) for x in pos])))
        for i,s in enumerate(samps):
            fout.write("%s\t%s\n" % (s, "\t".join([str(x) for x in z[i]])))


def plotInfo(info, samps):
    geneName, alPos, thisMin, thisMax = info
    covPos = range(thisMin, thisMax+1)
    posCounts = {n:{x:0 for x in covPos} for n in samps}

    for each in samps:
        enrPeps = io.fileList(each, header=False)
        for p in enrPeps:
            if p in alPos:
                for pos in alPos[p][2]:
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

