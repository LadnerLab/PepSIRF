#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:27:34 2021

@author: ANNAB
"""

#Import Modules
import pandas as pd
import optparse
import csv
import os
import inout as io    # https://github.com/jtladner/Modules/blob/main/inout.py

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def main():
    p = optparse.OptionParser()
    
    #Input controls
    p.add_option('-d', '--data',  help='Matrix containing data that will be used to generate correlations [none, REQ]')
    p.add_option('-s','--speciesD', default=False, help='Matirx containing data that will be used to find species names of peptides [none]')
    p.add_option('--nameHead', default="CodeName", help='Header name for column in metadata file that contains peptide names [CodeName]')
    p.add_option('--spHead', default="SpeciesID", help='Header name for column in metadata file that contains species info [SpeciesID]')

    
    #Correlation controls
    p.add_option('--minMax', type=float, default=30, help='Minimum maximum value [30]')
    p.add_option('-m', '--method', default='kendall', help='Method of correlation to compute pairwise correlation (kendall, pearson, spearman) [kendall]')
    p.add_option('--minCorr', type=float, default=0.5, help='Minimum correlation to be written in output file (0-1 inclusive) [0.5]')
    
    #Pairs controls
    p.add_option('--pairScatter', default=False, help='Generate scatter plots comparing the pairs of peptides above the threshold. Provide output directory [none]')
    
    #Output Controls
    p.add_option('-o', '--output', help='File name for output file [none, REQ]') 
    
    opts, args = p.parse_args()

    #Read in data file
    pyMat = pd.read_csv(opts.data, sep='\t', index_col=0)
    
    df = pd.DataFrame()
    
    #Get peptide names and values to add to dataframe
    for a,b in pyMat.iterrows():
        if max(b)>=opts.minMax:
            df[a] = b
    
    cmCorr = df.corr(method=opts.method)

    # Reformat matrix and remove any pairs that don't meet the correlation threshold
    dataCorr = cmCorr[abs(cmCorr) >= opts.minCorr].stack().reset_index()
    # Remove comparisons to self
    dataCorr = dataCorr[dataCorr['level_0'].astype(str)!=dataCorr['level_1'].astype(str)]
    dataCorr.rename(columns={"level_0":"peptide_A", "level_1":"peptide_B", 0:"%s-CorrCoef" % (opts.method)}, inplace=True)
    # Check to make sure the dataframe isn't empty
    if not dataCorr.empty:
        
        # filtering out lower/upper triangular duplicates 
        fn = lambda row: '-'.join(sorted([row['peptide_A'],row['peptide_B']]))
        col = dataCorr.apply(fn,axis=1)
        dataCorr = dataCorr.assign(orderedCols=col)
        dataCorr = dataCorr.drop_duplicates(['orderedCols'])
        dataCorr.drop('orderedCols', axis=1, inplace=True)

        # Add species info columns, if metadata provided
        if opts.speciesD:
            spD = io.fileDictHeader(opts.speciesD, opts.nameHead, opts.spHead)

            fnA = lambda row: spD[row['peptide_A']] 
            colA = dataCorr.apply(fnA,axis=1)
            dataCorr = dataCorr.assign(speciesA=colA)

            fnB = lambda row: spD[row['peptide_B']] 
            colB = dataCorr.apply(fnB,axis=1)
            dataCorr = dataCorr.assign(speciesB=colB)

    dataCorr.to_csv(opts.output, sep="\t", index=False)
    
    #Generate peptide pair scatters above given threshold
    if opts.pairScatter:
        
        #create dictionary of peptide data
        pD = df.to_dict('list')

        if not os.path.isdir(opts.pairScatter):
            os.mkdir(opts.pairScatter)
        else:
            print("Warning: The output directory %s already exists!" % opts.pairScatter)
        
        for index, row in dataCorr.iterrows():
        
            fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')
            
            pA = row["peptide_A"]
            pB = row["peptide_B"]

            x = pD[pA]
            y = pD[pB]

            ax.scatter(x, y, alpha=0.5, c='blue')
            ax.set_xlabel(pA, fontsize=15)
            ax.set_ylabel(pB, fontsize=15)
            fig.savefig('%s/%s~%s' % (opts.pairScatter, pA, pB), dpi=300, bbox_inches='tight')
            plt.close(fig)
                                            
    
                        
#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()