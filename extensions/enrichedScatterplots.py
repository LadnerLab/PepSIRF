#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import optparse, glob, os
import numpy as np
import inout as io             #Available at https://github.com/jtladner/Modules

# This script is used for generating scatterplots comparing control and experimental samples, with enriched peptides highlighted with a different color

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    #Experimental data
    p.add_option('-d', '--data',  help='Matrix containing data that will be used to generate scatterplots. [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")

    #Lists of enriched peptides
    p.add_option('-e', '--enrDir',  help='Directory containing lists of enriched peptides. [None, REQ]')
    p.add_option('-f', '--enrFile', help='Tab delimited file containing data that will be used to generate scatterplots. [None, OPT]')
    p.add_option('-z', '--zenrFiles', help='Tab delimited file containing z score thresholds and related directories containing data that will be used to generate scatterplots. [None, OPT]')
    p.add_option('--enrExt', help='Common file ending for files containing enriched sets of peptides. Once removed, the remaining filename should consist only of sample name(s) [None, REQ]')
    p.add_option('--snDelim', default="~", help='Delimiter used to separate sample names in pEnrich files [~]')

    #Negative contols
    p.add_option('--negMatrix',  help='Optional way to provide a separate data matrix for negative controls. If not provided, negative controls will be assumed to be from the same matrix as the experiemntal data  [None, OPT]')
    p.add_option('-c', '--negControls',  help='Comma-delimited names of samples to use as negative controls. If this is not provided, then all of the samples in the negative control matrix will be used [None, REQ]')
    p.add_option('-i', '--negative_id', help='Optional approach for identifying negative controls. Provide a unique string at the start of all negative control samples. [None, OPT]')
    p.add_option('-x', '--xLabel', default="", help='Label to be used for x-axis. []')

    #Output controls
    p.add_option('-o', '--outDir', help='Directory name for output files. Will be created. [None]')
    p.add_option('--plotType', default="png", help='Type of plot to generate. This should be a file extension recognized by matplotlib, e.g., pdf, png, tiff [png]')
    p.add_option('--plotLog', type=float, default=False, help="Use if you want axes to be shown on a log-scale. Argument provided should be a float to add to the axes before calcluating the log value [False]")
    p.add_option('--negColor', default="#1b9e77", help='Color to use for Unenriched peptides. [#1b9e77]')
    p.add_option('--posColor', default="#d95f02", help='Color to use for Enenriched peptides. [#d95f02]')
    p.add_option('--ptsSize', type=float, default=10, help='Size to use for plotted points. [10]')
    p.add_option('--ptsTrans', type=float, default=0.5, help='Transparency to use for plotted points. [0.5]')

    opts, args = p.parse_args()
    
    #Check if enriched directory and file were both provided and issue warning
    if opts.enrFile and opts.enrDir:
        print("Warning: An enriched directory and an enriched file were provided. File will be used, and directory will be ignored.")

    #Create output dir
    if os.path.isdir(opts.outDir):
        print("The output diretory %s already exists!" % (opts.outDir))
    else:
        os.mkdir(opts.outDir)

    #Read in data file 
    dataD = io.fileDictFullRowNames(opts.data, opts.delim)

    #Prep negative control data for plotting on x-axis, this will be the same across all of the plots
    if opts.negMatrix:
        negD = io.fileDictFullRowNames(opts.negMatrix, opts.delim)
    else:
        negD = dataD
    
    if opts.negative_id:
        negNames = []
        for name in negD.keys():
            if name.startswith(opts.negative_id):
                negNames.append(name)
        if not negNames:
            print("Warning: No samples found matching the provided negative ID")
    elif opts.negControls:
        negNames = opts.negControls.split(",")
    else:
        if opts.negMatrix:
            print("Warning: Using all samples in %s to generate x-axis values" % opts.negMatrix)
        else:
            print("Warning: Using all samples in %s to generate x-axis values" % opts.data)
        negNames = list(negD.keys())
    
    peptideNames = list(negD[negNames[0]].keys())
    
    x = [np.mean([float(negD[sn][pn]) for sn in negNames]) for pn in peptideNames]
    if opts.plotLog:
        x = [np.log10(point+opts.plotLog)for point in x]
    
    # Read in lists of enriched peptides
    if opts.enrFile:
        filePath = {}
        with open(opts.enrFile, 'r') as fin:
            for line in fin:
                filelist = line.strip("\n").split("\t")
                filePath[filelist[0]] = []
                i = 1
                while i < len(filelist):
                    if filelist[i] =='':
                        i += 1
                    else:
                        filePath[filelist[0]].append(filelist[i])
                        i += 1
        enrFiles = list(filePath.keys())
    elif opts.enrDir:
        enrFiles = glob.glob("%s/*%s" % (opts.enrDir, opts.enrExt))

    #Generate plot for each list of enriched peptides
    if opts.enrFile or opts.enrDir:
        for eF in enrFiles:
            enrichedD = io.fileEmptyDict(eF, header=False)
            if opts.enrFile:
                sNames = filePath[eF]
            elif opts.enrDir:
                sNames = os.path.basename(eF).split(opts.enrExt)[0].split(opts.snDelim)
            
            # Remove any sample names that are not found in the dataD
            sNames = [s for s in sNames if s in dataD]
            
            if len(sNames) > 0:
            
                #Average values for y-axis
                y = [np.mean([float(dataD[sn][pn]) for sn in sNames]) for pn in peptideNames]
            
                #Determine color for each point
                c = [opts.posColor if p in enrichedD else opts.negColor for p in peptideNames]
            
                #Generate plot
                fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w') 
            
                if opts.plotLog:
                    y = [np.log10(point+opts.plotLog)for point in y]
                    ax.set_ylabel(",".join(sNames)+f" log10(value+{opts.plotLog})", fontsize=15)
                    ax.set_xlabel(opts.xLabel+f" log10(value+{opts.plotLog})", fontsize=15)
                else:
                    ax.set_ylabel(",".join(sNames), fontsize=15)
                    ax.set_xlabel(opts.xLabel, fontsize=15)  
                
                ax.scatter(x, y, c=c, alpha=opts.ptsTrans, s=opts.ptsSize)
    
        #        if lim:
        #            ax.set_xlim(lim[0], lim[1])
        #            ax.set_ylim(lim[0], lim[1])
    
                fig.savefig("%s/%s.%s" % (opts.outDir, opts.snDelim.join(sNames), opts.plotType), dpi=300, bbox_inches='tight')
    
                plt.close(fig)
            
            else:
                print("Skipping %s because none of the associated replicates were found in %s" % (eF, opts.data))
                
    elif opts.zenrFiles:
        #create lists to collect z score thresholds and directory paths
        zThresh = {}
        
        #create list for colors of different z score thresholds
        clrs = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000']
        
        #open file to collect data from input file
        with open(opts.zenrFiles, 'r') as fin:
            for line in fin:
                zT,dP = line.strip("\n").split("\t")
                zThresh[int(zT)] = []
                zThresh[int(zT)].append(dP)      
                
        #collect z score thresholds in ascending order
        zScore = list((sorted(zThresh.keys())))

        #collect enriched file names
        enrichedD = {}
        for zS in zThresh:
            for dP in zThresh[zS]:
                enrFiles = glob.glob("%s/*%s" % (dP, opts.enrExt))
                for eF in enrFiles:
                    enrichedD[os.path.basename(eF)] = []

        #collect file paths to the enriched files for each z score theshold in ascending order
        for zS in zScore:
            for dP in zThresh[zS]:
                enrFiles = glob.glob("%s/*%s" % (dP, opts.enrExt))
                for eF in enrFiles:
                    for file in enrichedD:
                        if os.path.basename(eF) == file:
                            enrichedD[os.path.basename(eF)].append(eF)
                            
        #generate plots for each enriched file
        for file in enrichedD:
            
            sNames = file.split(opts.enrExt)[0].split(opts.snDelim)
                
            # Remove any sample names that are not found in the dataD
            sNames = [s for s in sNames if s in dataD]
            
            if len(sNames) > 0:
            
                #Average values for y-axis
                y = [np.mean([float(dataD[sn][pn]) for sn in sNames]) for pn in peptideNames]
                
                if opts.plotLog:
                        y = [np.log10(point+opts.plotLog)for point in y]
            
            #generate plots
            fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')

            posC = 0
            c = ['#808080' for p in peptideNames]
            
            #add legend title
            ax.text(0.87,0.95,'Z Threshold',fontsize=10,ha='center', va='center_baseline', transform=ax.transAxes)
            
            #parse through every enriched file for each z score theshold
            for path in enrichedD[file]:
                peptideD = io.fileEmptyDict(path, header=False)
                    
                #Determine color for each point
                i = 0
                for p in peptideNames:
                    if p in peptideD:
                        c[i] = clrs[posC]
                        i+=1
                    else:
                        i+=1
                posC += 1

            #add axis labels
            if opts.plotLog:
                ax.set_ylabel(",".join(sNames)+f" log10(value+{opts.plotLog})", fontsize=15)
                ax.set_xlabel(opts.xLabel+f" log10(value+{opts.plotLog})", fontsize=15)
            else:
                ax.set_ylabel(",".join(sNames), fontsize=15)
                ax.set_xlabel(opts.xLabel, fontsize=15)
            
            #add legend items in descending order
            txtC = posC - 1
            txtA = 0.90
            for zS in list(reversed(zScore)):
                ax.text(0.87,txtA,zS,c=clrs[txtC],fontsize=10,ha='center', va='center_baseline', transform=ax.transAxes)
                txtC-=1
                txtA-=0.05
                
            #plot points and save figure
            ax.scatter(x, y, c=c, alpha=opts.ptsTrans, s=opts.ptsSize)
            fig.savefig("%s/%s.%s" % (opts.outDir, opts.snDelim.join(sNames), opts.plotType), dpi=300, bbox_inches='tight')
            plt.close(fig)
    
    else:
        print("Warning: No enriched directory, file, or z file provided. Must provide one.")
#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

