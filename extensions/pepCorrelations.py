# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:27:34 2021

@author: ANNAB
"""

#Import Modules
import pandas as pd
import optparse
import csv

def main():
    p = optparse.OptionParser()
    
    #Input controls
    p.add_option('-d', '--data',  help='Matrix containing data that will be used to generate correlations [none, REQ]')
    
    #Correlation controls
    p.add_option('--minMax', type=float, default=100, help='minimum maximum value [100]')
    p.add_option('-m', '--method', default='kendall', help='method of correlation to compute pairwise correlation (kendall, pearson, spearman) [kendall]')
    p.add_option('--minCorr', type=float, default=0.5, help='Minimum correlation to be written in output file (0-1 inclusive) [0.5]')
    
    #Output Controls
    p.add_option('-o', '--output', help='File name for output file [none, REQ]') 
    p.add_option('-s','--speciesD', default=False, help='Matirx containing data that will be used to find species names of peptides [none]')
    
    opts, args = p.parse_args()
    
    #Read in data file
    pyMat = pd.read_csv(opts.data, sep='\t', index_col=0)
    
    df = pd.DataFrame()
    
    #Get peptide names and values to add to dataframe
    for a,b in pyMat.iterrows():
        if max(b)>=opts.minMax:
            df[a] = b
    
    cmCorr = df.corr(method=opts.method)
    
    if opts.speciesD:
        #Read in species data file
        pyMat = pd.read_csv(opts.speciesD, sep='\t', index_col=0, low_memory=False)
    
        df = pd.DataFrame()
    
        #get peptide information
        for pep1,corrS in cmCorr.iteritems():
            for a,b in pyMat.iterrows():
                if a == pep1:
                    df[a] = b

    #Generate tab-delimited correlation file
    with open(opts.output, 'w', newline='') as test_out:
        tsv_writer = csv.writer(test_out, delimiter='\t')
        
        if not opts.speciesD:
            tsv_writer.writerow(['peptide1', 'peptide2', 'correlation('+opts.method+')'])
        else:
            tsv_writer.writerow(['peptide1', 'peptide2', 'species1', 'species2', 'correlation('+opts.method+')'])
        
        for pep1,corrS in cmCorr.iteritems():
            for pep2, corr in corrS.iteritems():
                #Minimum correlation to be added to output file
                if corr >= opts.minCorr:
                    
                    #Write row to output file, excluding comparisons against self
                    if pep1 != pep2:
                        if not opts.speciesD:
                            tsv_writer.writerow([pep1, pep2, corr])
                        
                        #write rows containing species names to output file
                        else:
                            for pep in df:
                                if pep1 == pep:
                                    species1 = df[pep]['Species']
                                    
                                if pep2 == pep:
                                    species2 = df[pep]['Species']
                                    
                            tsv_writer.writerow([pep1, pep2, species1, species2, corr])
                        
#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()