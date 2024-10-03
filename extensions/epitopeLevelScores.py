#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import argparse, os
import pandas as pd
pd.options.mode.copy_on_write = True
import seaborn as sns
import numpy as np
import inout as io


# This script reads in 1) normalized or raw read counts and 2) bins and calculated zscores

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-d', '--data',  help='Data matrix containing scores of interest.', required=True)

    parser.add_argument('--outDir', default=".",  help='Where to generate output files.')
    parser.add_argument('--delim', default="_", help="Delimiter used in the sample names.")
    parser.add_argument('--catIndex', default=1, type=int, help="Index to use to extract catgeory from sample name, after splitting on --delim.")
#    parser.add_argument('--sumStatOut', help="Name for file that will contain multi-epitope summary statistics.")
    parser.add_argument('--sumStatNeg', help="Name for negative control category to use in calculating summary statistics.")
    parser.add_argument('rest', nargs=argparse.REMAINDER)
    
    opts = parser.parse_args()

    # Read in data file
    sD = io.fileDictFullRowNames(opts.data, valType="float")
    
    epiDir = "/Users/jtladner/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/LL_Projects/RAPTER/RPTR1/mappingEpitopes"

    # Generate max Z score output
    for epiPepF in opts.rest:
        epiPepD = io.fileDictHeader(epiPepF, "Epitope", "Peptides")
        epiPepD = {k:v.split("~") for k,v in epiPepD.items()}
        
        sOut = f"{opts.outDir}/epiMax_{os.path.basename(epiPepF)}"
        with open(sOut, "w") as fout:
            fout.write("Sample\tCategory\tEpitope\tMaxZ\n")
            for en, epL in epiPepD.items():
                for sn, thisZ in sD.items():
                    cn = sn.split(opts.delim)[opts.catIndex]
                    mz = max([thisZ[pn] for pn in epL])
                    fout.write(f"{sn}\t{cn}\t{en}\t{mz}\n")

        # Generate epitope level maxZ boxplot
        df = pd.read_csv(sOut, sep="\t", header=0)
    
        # Generate boxplot figure
        fig,ax = plt.subplots(1, 1, figsize=(len(epiPepD),5),facecolor='w')
        sns.boxplot(data=df, x="Epitope", hue="Category", y="MaxZ", ax=ax, fliersize=0)
        sns.stripplot(data=df, x="Epitope", hue="Category", y="MaxZ", jitter=True, dodge=True, linewidth=0.5, ax=ax, zorder=3)
        
        #Modify axis labels
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_ha('right')

        
        #Modify legend
        handles, labels = ax.get_legend_handles_labels()
        numCats = int(len(handles)/2)
        plt.legend(handles[0:numCats], labels[0:numCats], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig.savefig(f'{sOut}.png',dpi=300,bbox_inches='tight')
        fig.savefig(f'{sOut}.pdf',dpi=300,bbox_inches='tight')


        # Generate epitope level maxZ boxplot, log2 transformed
        df = pd.read_csv(sOut, sep="\t", header=0)
    
        # Generate boxplot figure
        fig,ax = plt.subplots(1, 1, figsize=(len(epiPepD),5),facecolor='w')
        sns.boxplot(x=df["Epitope"], hue=df["Category"], y=df["MaxZ"]+1, ax=ax, fliersize=0, log_scale=2)
        sns.stripplot(x=df["Epitope"], hue=df["Category"], y=df["MaxZ"]+1, jitter=True, dodge=True, linewidth=0.5, ax=ax, zorder=3, log_scale=2)
        
        #Modify axis labels
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_ha('right')

        
        #Modify legend
        handles, labels = ax.get_legend_handles_labels()
        numCats = int(len(handles)/2)
        plt.legend(handles[0:numCats], labels[0:numCats], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig.savefig(f'{sOut}_log2.png',dpi=300,bbox_inches='tight')
        fig.savefig(f'{sOut}_log2.pdf',dpi=300,bbox_inches='tight')

        

        # Generate multi-peptide summary stats, if requested
        if opts.sumStatNeg:
            # Set name of summary stat output file
            opts.sumStatOut = f"epiMaxSummaryStats_{os.path.basename(epiPepF)}"
        
            s2catD = {r.Sample:r.Category for i,r in df.iterrows()} 
            predf = {"Sample":[], "Category":[]}
            for k,v in s2catD.items():
                predf["Sample"].append(k)
                predf["Category"].append(v)
            odf = pd.DataFrame(predf, index=predf["Sample"])
            odf["Enriched Epitopes"] = 0
            odf["Zscore Sum"] = 0.0

            for epi in set(df["Epitope"]):
                epiDF = df[df["Epitope"]==epi]
    
                #Pull out controls and calculate threshold
                controls = epiDF[epiDF["Category"]==opts.sumStatNeg]
                zThresh = np.mean(controls["MaxZ"]) +  2*np.std(controls["MaxZ"])
    
                #Step through all samples and collect info
                for i,row in df.iterrows():
                    if row.Epitope==epi:
                        if row.MaxZ >zThresh:
                            odf.loc[row.Sample, "Enriched Epitopes"]+=1
                            odf.loc[row.Sample, "Zscore Sum"]+=row.MaxZ
            
                # Write out summary stat file
                odf.to_csv(opts.sumStatOut, sep="\t", index=False)
            
                    
            # Generate boxplot figure
            fig,ax = plt.subplots(1, 2, figsize=(12,5),facecolor='w')
            
            # Z score sum panel
            sns.boxplot(data=odf, hue="Category", y="Zscore Sum", ax=ax[0], fliersize=0)
            sns.stripplot(data=odf, hue="Category", y="Zscore Sum", jitter=True, dodge=True, linewidth=0.5, ax=ax[0], zorder=3, legend=False)

            # Enriched epitopes panel
            sns.boxplot(data=odf, hue="Category", y="Enriched Epitopes", ax=ax[1], fliersize=0)
            sns.stripplot(data=odf, hue="Category", y="Enriched Epitopes", jitter=True, dodge=True, linewidth=0.5, ax=ax[1], zorder=3, legend=False)

            for label in ax[0].get_xticklabels() + ax[1].get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')

            # Move legends
            sns.move_legend(ax[0], "lower center", bbox_to_anchor=(0.5, -0.35))
            sns.move_legend(ax[1], "lower center", bbox_to_anchor=(0.5, -0.35))

            #Save figure
            fig.savefig(f'{opts.sumStatOut}.png',dpi=300,bbox_inches='tight')
            fig.savefig(f'{opts.sumStatOut}.pdf',dpi=300,bbox_inches='tight')


            # Generate log2 boxplot figure
            fig,ax = plt.subplots(1, 2, figsize=(12,5),facecolor='w')
            
            # Z score sum panel
            sns.boxplot(hue=odf["Category"], y=odf["Zscore Sum"]+1, ax=ax[0], fliersize=0, log_scale=2)
            sns.stripplot(hue=odf["Category"], y=odf["Zscore Sum"]+1, jitter=True, dodge=True, linewidth=0.5, ax=ax[0], zorder=3, log_scale=2, legend=False)

            # Enriched epitopes panel
            sns.boxplot(data=odf, hue="Category", y="Enriched Epitopes", ax=ax[1], fliersize=0)
            sns.stripplot(data=odf, hue="Category", y="Enriched Epitopes", jitter=True, dodge=True, linewidth=0.5, ax=ax[1], zorder=3, legend=False)

            for label in ax[0].get_xticklabels() + ax[1].get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')
            
            # Move legends
            sns.move_legend(ax[0], "lower center", bbox_to_anchor=(0.5, -0.35))
            sns.move_legend(ax[1], "lower center", bbox_to_anchor=(0.5, -0.35))
            
            #Save figure
            fig.savefig(f'{opts.sumStatOut}_log2.png',dpi=300,bbox_inches='tight')
            fig.savefig(f'{opts.sumStatOut}_log2.pdf',dpi=300,bbox_inches='tight')


#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

