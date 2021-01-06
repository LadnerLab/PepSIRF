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
    p.add_option('--enrExt', help='Common file ending for files containing enriched sets of peptides. Once removed, the remaining filename should consist only of sample name(s) [None, REQ]')
    p.add_option('--snDelim', default="~", help='Delimiter used to separate sample names in pEnrich files [~]')

    #Negative contols
    p.add_option('--negMatrix',  help='Optional way to provide a separate data matrix for negative controls. If not provided, negative controls will be assumed to be from the same matrix as the experiemntal data  [None, OPT]')
    p.add_option('-c', '--negControls',  help='Comma-delimited names of samples to use as negative controls. If this is not provided, then all of the samples in the negative control matrix will be used [None, REQ]')
    p.add_option('-x', '--xLabel', default="", help='Label to be used for x-axis. []')

    #Output controls
    p.add_option('-o', '--outDir', help='Directory name for output files. Will be created. [None]')
    p.add_option('--plotType', default="png", help='Type of plot to generate. This should be a file extension recognized by matplotlib, e.g., pdf, png, tiff [png]')
    p.add_option('--plotLog', default=False, action="store_true", help="Use if you want axes to be shown on a log-scale. [False]")
    p.add_option('--negColor', default="#1b9e77", help='Color to use for Unenriched peptides. [#1b9e77]')
    p.add_option('--posColor', default="#d95f02", help='Color to use for Enenriched peptides. [#d95f02]')

    opts, args = p.parse_args()

    #Create output dir
    if os.path.isdir(opts.outDir):
        print("The output direory %s already exists!" % (opts.outDir))
    else:
        os.mkdir(opts.outDir)

    #Read in data file 
    dataD = io.fileDictFullRowNames(opts.data, opts.delim)
    
    #Prep negative control data for plotting on x-axis, this will be the same across all of the plots
    if opts.negMatrix:
        negD = io.fileDictFullRowNames(opts.negMatrix, opts.delim)
    else:
        negD = dataD
    
    if not opts.negControls:
        negNames = list(negD.keys())
    else:
        negNames = opts.negControls.split(",")
    
    peptideNames = list(negD[negNames[0]].keys())
    
    x = [np.mean([float(negD[sn][pn]) for sn in negNames]) for pn in peptideNames]

    
    # Read in lists of enriched peptides
    enrFiles = glob.glob("%s/*%s" % (opts.enrDir, opts.enrExt))

    #Generate plot for each list of enriched peptides
    for eF in enrFiles:
        enrichedD = io.fileEmptyDict(eF, header=False)
        sNames = eF.split("/")[-1].split(opts.enrExt)[0].split(opts.snDelim)
        
        #Average values for y-axis
        y = [np.mean([float(dataD[sn][pn]) for sn in sNames]) for pn in peptideNames]
        
        #Determine color for each point
        c = [opts.posColor if p in enrichedD else opts.negColor for p in peptideNames]
        
        #Generate plot
        fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')            
        
        if opts.plotLog:
            ax.scatter(np.log10(x), np.log10(y), c=c, alpha=0.5)
        else:
            ax.scatter(x, y, c=c, alpha=0.5)

        ax.set_ylabel(",".join(sNames), fontsize=15)
        ax.set_xlabel(opts.xLabel, fontsize=15)

#        if lim:
#            ax.set_xlim(lim[0], lim[1])
#            ax.set_ylim(lim[0], lim[1])

        fig.savefig("%s/%s.%s" % (opts.outDir, opts.snDelim.join(sNames), opts.plotType), dpi=300, bbox_inches='tight')


#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

