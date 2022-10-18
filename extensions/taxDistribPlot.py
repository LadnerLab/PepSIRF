#!/usr/bin/env python

import argparse, glob, os
import inout as io
import seaborn as sns

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.patheffects as path_effects
from IPython.display import HTML
from matplotlib.patches import Rectangle

typeface='Arial'
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.family']=typeface
mpl.rcParams['font.size']=22
mpl.rcParams['hatch.linewidth'] = 0.4


from collections import defaultdict


# The purpose of this script is to generate a figure summarizing hits across several related taxomonic groups

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--m1_pepName', default="CodeName", help="Header linked to peptide name column in metadata file 1.")
    parser.add_argument('--m1_sid', default="SpeciesID", help="Header linked to species ID column in metadata file 1.")
    parser.add_argument('--m1_cat', default="Category", help="Header linked to Category column in metadata file 1.")
    parser.add_argument('--m2_sid', default="Species ID", help="Header linked to species ID column in metadata file 2.")
    parser.add_argument('--m2_sp', default="Species", help="Header linked to species name column in metadata file 2.")
    parser.add_argument('--m2_gen', default="Genus", help="Header linked to genus name column in metadata file 2.")
    parser.add_argument('--m2_fam', default="Family", help="Header linked to family name column in metadata file 2.")
    parser.add_argument('--enrFileStr', default="*_enriched.txt", help="String with wildcards that will be used to glob enriched peptide files. This should NOT include file paths.")
    parser.add_argument('--cats2incl', default="SetCover", help="Comma-separated list of categories to include in counts.")
    parser.add_argument('--fnSep', default="_", help="Separator used to divide pieces of info from enriched peptide filenames.")
    parser.add_argument('--sampNamePos', default=0, type=int, help="Position in enriched peptide filenames that corresponds to ther sample name.")
    parser.add_argument('--fontSize1', default=20, type=int, help="Font size for primary labels.")
    parser.add_argument('--fontSize2', default=15, type=int, help="Font size for secondary labels.")
    parser.add_argument('--cmap', default="rainbow", help="Name of colormap to use for plot.")


    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('--m1', help='Metadata file used to link peptide names to species IDs and categories', required=True)
    reqArgs.add_argument('--m2', help='Metadata file used to species IDs to species, genus and family names', required=True)
    reqArgs.add_argument('-d', '--dir',  help='Directory containing files with lists of enriched peptides.', required=True)
    reqArgs.add_argument('-o', '--out',  help='Name for output graph, including desired file extension.', required=True)
    reqArgs.add_argument('--t1',  help='Broadest taxonomic category for inclusion in the plot. Normally of lower resolution than t2.', choices=["species", "genus", "family"], required=True)
    reqArgs.add_argument('--t2',  help='Taxonomic category that will actually be plotted. Normally of higher resolution than t1.', choices=["species", "genus", "family"], required=True)
    reqArgs.add_argument('--n1',  help='Taxonimic name corresponding to t1.', required=True)

    args = parser.parse_args()
    
    #Parse categories of interest
    cats2incl = args.cats2incl.split(",")
    
    # Read in required metadata
    sidD = io.fileDictHeader(args.m1, args.m1_pepName, args.m1_sid)
    catD = io.fileDictHeader(args.m1, args.m1_pepName, args.m1_cat)

    spD = io.fileDictHeader(args.m2, args.m2_sid, args.m2_sp)
    genD = io.fileDictHeader(args.m2, args.m2_sid, args.m2_gen)
    famD = io.fileDictHeader(args.m2, args.m2_sid, args.m2_fam)
    
    # Parse lists of enriched peptides
    enrFL = glob.glob("%s/%s" % (args.dir, args.enrFileStr))
    enrPepLists = defaultdict(list)

    # Make lists of enriched peptides for each sample, pooled across capture proteins, when relevant
    for each in enrFL:
        samp = os.path.basename(each).split(args.fnSep)[args.sampNamePos]
        enrPepLists[samp] += [x for x in io.fileList(each, header=False) if catD[x] in cats2incl]
    
    # Convert to set and then back to list to prevent double counting of the same peptide from different capture proteins
    for smp, pList in enrPepLists.items():
        enrPepLists[smp] = list(set(pList))
    
    # Generate counts of the number of peptides per category of interest
    totalPepsPer = defaultdict(int)

    for cn, cat in catD.items():
        if cat in cats2incl:
            sid = sidD[cn]

            if sid in genD:

                taxNames = {
                    "species": spD[sid],
                    "genus": genD[sid],
                    "family": famD[sid],
                    }
                    
                if taxNames[args.t1] == args.n1:
                    focal = taxNames[args.t2]
                    totalPepsPer[focal]+=1



    # Count number of peptides from each focal group in each sample
    tPepLD = defaultdict(list)

    for samp, pepL in enrPepLists.items():
        thisD = defaultdict(int)
        for p in pepL:
            sid = sidD[p]
            if sid in genD:

                taxNames = {
                    "species": spD[sid],
                    "genus": genD[sid],
                    "family": famD[sid],
                    }

                if taxNames[args.t1] == args.n1:
                    focal = taxNames[args.t2]
                    thisD[focal]+=1

        for k,v in thisD.items():
            tPepLD[k].append(v)


    #Generate order for groups in figure and store number total peptides
    pepsInDesign = {}
    sl = sorted([(v,k) for k,v in totalPepsPer.items()])[::-1]
    for v,k in sl:
        pepsInDesign[k]=v
    
    #Generate colors to use in plot
    cmap = mpl.cm.get_cmap(args.cmap)
    colorD={}
    if len(pepsInDesign) == 1:
        increment = 1
    else:
        increment = 1/(len(pepsInDesign)-1)
    for i,k in enumerate(list(pepsInDesign.keys())):
        colorD[k] = cmap(0+i*increment)
    
    #Generate Figure
    fig,ax = plt.subplots(1,2,figsize=(4+2*len(totalPepsPer),4),facecolor='w',gridspec_kw={'width_ratios': [8, 1]})

    w=1
    x=[]
    y=[]
    bot=0

    for f in pepsInDesign:
        #Generating bar graphs
        ax[1].bar(args.n1, pepsInDesign[f], w, bottom=bot, color=colorD[f], edgecolor="k", linewidth=0.5)
        bot+=pepsInDesign[f]
        for val in tPepLD[f]:
            x.append("%s\n(%d)" % (f, len(tPepLD[f])))
            y.append(val)


    sns.boxplot(x=x, y=y, ax=ax[0], fliersize=0, palette={"%s\n(%d)" % (f, len(tPepLD[f])):c for f,c in colorD.items()})
    sns.stripplot(x=x, y=y, ax=ax[0], jitter=True, dodge=True, linewidth=0.5, palette={"%s\n(%d)" % (f, len(tPepLD[f])):c for f,c in colorD.items()})

    # Formatting boxplots
    ax[0].set_ylabel("Enriched Peptides", fontsize=args.fontSize1)
    ax[0].tick_params(axis='x', labelsize=args.fontSize2)

    # Formatting for bar graphs
    ax[1].set_ylim(0,bot)
    ax[1].set_xlim(-0.45,0.45)
    ax[1].set_ylabel("Assay Peptides", fontsize=args.fontSize1, labelpad=10)
    ax[1].tick_params(axis='x', labelsize=args.fontSize2)
    ax[1].yaxis.set_label_position("right")
    ax[1].yaxis.tick_right()

    plt.savefig(args.out,dpi=200,bbox_inches='tight')

#----------------------End of main()

    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

