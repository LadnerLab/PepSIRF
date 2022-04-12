#!/usr/bin/env python

import argparse, re, subprocess, glob, os
import numpy as np
import inout as io
from collections import defaultdict
import itertools as it

# If matplotlib is available, generate some summary figures
try: 
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt
    matplotReady=True
        
    
except:
    print("Not generating figures because of error. Check to see if matplotlib is available.")
    matplotReady=False

# Used to automate several consecutive analyses with PepSIRF and to summarize the results

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-e", "--binary", help="How to call the pepSIRF binary/executable", default="pepsirf")

    inArgs = parser.add_argument_group('input files')
    inArgs.add_argument("-c", "--colsum", help="Input colsum normalized count matrix.", required=True)
    inArgs.add_argument("-z", "--zscore", help="Z score matrix generated using buffer-only negative control bins.", required=True)
    inArgs.add_argument("-p", "--pairs", help="Tab-delimited file containing longitudinal sample pairs. The 1st column should be the name of the 1st timepoint, the 2nd column the name of the 2nd timepoint. Replicates can be provided as comma-separated lists.", required=True)
    inArgs.add_argument("-r", "--raw", help="Input raw count matrix. Required if you want to use a raw count threshold for enrich.")
#    inArgs.add_argument("-t", "--thresh", help="Threshold file for enrich. This is an optional starting place.")
#    inArgs.add_argument("-d", "--diff", help="Input diff normalized count matrix. This is an optional starting place.")
#    inArgs.add_argument("-f", "--diffratio", help="Input diff_ratio normalized count matrix. This is an optional starting place.")
#    inArgs.add_argument("-b", "--bins", help="Peptide bin file. Required for calculating Z scores.")
#    inArgs.add_argument("-n", "--names", help="File containing sample names. This is an optional starting place.")

    
    threshArgs = parser.add_argument_group('thresholds')
    threshArgs.add_argument("--zThresh_Neg", default="10,10", help="Threshold for Z scores generated relative to buffer-only negative controls. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument("--zThresh_t1", default="10,10", help="Threshold for Z scores generated relative to timepoint 1. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument("--csThresh", default="0,0", help="Colum-Sum norm count threshold. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument('--rawThresh',  default="15000", help='Total raw read count for a sample to be included in enrichment analyses. Can include up to two floating points separated by a comma.')
    threshArgs.add_argument('--hdi',  default=0.98, type=float, help='The highest density interval to be used for calculation of mean and stdev in the zscore module.')
#    threshArgs.add_argument("--sbdrThresh",  default=False, help="Negative control ratio threshold. Can include up to two floating points separated by a comma.")

#    controlArgs = parser.add_argument_group('control info')
#    controlArgs.add_argument("--negNormMat", help="Alternative colsum normalized matrix form which to obtain the negative controls.")
#    controlArgs.add_argument("--negative_id", help="Optional approach for identifying negative controls. Provide a unique string at the start of all negative control samples.")
#    controlArgs.add_argument("--negative_names", help="Optional approach for identifying negative controls. Comma-separated list of negative control sample names.")

#     enrichArgs = parser.add_argument_group('enrich options')
#     enrichArgs.add_argument("--v13", default=False, action="store_true", help="Flag to be used if using pepsirf version 1.3 to run p_ernich module. If not set, pepsirf version 1.4 will be assumed, and enrich module will be run.")
#     enrichArgs.add_argument("--sEnrich", default=False, action="store_true", help="Generate lists of enriched peptides separately for each pulldown. Will actually run p_enrich/enrich, but with the same sample specified for each replicate.")
#     enrichArgs.add_argument("--inferPairs", default=False, action="store_true", help="Infer sample pairs from names. This option assumes names of replicates will be identical with the exception of a final string denoted with a '_'. For example, these names would be considered two replicates of the same sample: VW_100_1X_A and VW_100_1X_B")

    outArgs = parser.add_argument_group('output options')
    outArgs.add_argument("-o", "--outDir", help="Name for base output directory.", required=True)
#     outArgs.add_argument("--repZscatters", default=False, action="store_true", help="Generate scatter plots comparing Z scores for all sample pairs. Even if set at command line, will be turned off when '--sEnrich' is used.")
#     outArgs.add_argument("--repCSscatters", default=False, action="store_true", help="Generate scatter plots comparing colu_sum normalized scores for all sample pairs. Even if set at command line, will be turned off when '--sEnrich' is used.")
#     outArgs.add_argument("--scatterFormat", default="png", help="Output file format for replicate scatterplots.")
#     outArgs.add_argument("--enrichedScatters", default=False, help="Generate scatter plots comparing col-sum normalized read counts between a sample and negative controls. Argument provided should mirror the way that should enrichedScatterplots.py  be called on your system.")

    args = parser.parse_args()
    
    # Read in sample pairs to analyze
    pairD = io.fileDict(args.pairs, key=0, val=1, delim="\t", header=False)
    pairD = {tuple(k.split(",")):tuple(v.split(",")) for k,v in pairD.items()}
    
    # Make directories that will be needed for the output files
    io.makeDir(args.outDir)
    io.makeDir("%s/nameLists" % args.outDir)
    io.makeDir("%s/nameLists/repFiles" % args.outDir)
    io.makeDir("%s/t1_CS" % args.outDir)
    io.makeDir("%s/pair_CS" % args.outDir)
    io.makeDir("%s/bins" % args.outDir)
    io.makeDir("%s/pair_diff" % args.outDir)
    io.makeDir("%s/pair_Z" % args.outDir)
    io.makeDir("%s/threshFiles" % args.outDir)
#    io.makeDir("%s/enrich" % args.outDir)
    io.makeDir("%s/deconv" % args.outDir)
    
    # Step through each pair of samples to process
    for t1, t2 in pairD.items():
        base = "%s_%s" % (sharedSubStr(t1), sharedSubStr(t2))
        
        # Generate name lists that will be used by subjoin
        t1_nameF = "%s/nameLists/%s_t1.txt" % (args.outDir, base)
        t1t2_nameF = "%s/nameLists/%s.txt" % (args.outDir, base)
        io.writeList(t1, t1_nameF, delim="\n")
        io.writeList(t1+t2, t1t2_nameF, delim="\n")
        
        # Generate subset colsum files
        t1_subCS = "%s/t1_CS/%s_t1_CS.tsv" % (args.outDir, base)
        cmd = '%s subjoin -i "%s",%s -o %s >> norm.out' % (args.binary, args.colsum, t1_nameF, t1_subCS)
        print(cmd)
        subprocess.run(cmd, shell=True)
        
        t1t2_subCS = "%s/pair_CS/%s_CS.tsv" % (args.outDir, base)
        cmd = '%s subjoin -i "%s",%s -o %s >> norm.out' % (args.binary, args.colsum, t1t2_nameF, t1t2_subCS)
        print(cmd)
        subprocess.run(cmd, shell=True)
        
        # Generate sample/timepoint specific bins
        binF = "%s/bins/%s_min300r1_bins.tsv" % (args.outDir, base)
        cmd = "%s bin -s %s -b 300 -r 1 -o %s >> bin.out" % (args.binary, t1_subCS, binF)
        print(cmd)
        subprocess.run(cmd, shell=True)

        # Generate sample/timepoint specific diff measures
        diffF = "%s/pair_diff/%s_Diff.tsv" % (args.outDir, base)
        cmd = '%s norm -a diff -p %s -o %s --negative_control "%s" --negative_names %s >> norm.out' % (args.binary, t1t2_subCS, diffF, t1_subCS, ",".join(t1))
        print(cmd)
        subprocess.run(cmd, shell=True)

        # Generate timepoint specific Z scores
        zF = "%s/pair_Z/%s_Z-HDI%d.tsv" % (args.outDir, base, int(args.hdi*100))
        zNan = "%s/pair_Z/%s_Z-HDI%d.nan" % (args.outDir, base, int(args.hdi*100))
        cmd = '%s zscore -s %s -o %s -n %s -b "%s" -d %f >> zscore.out' % (args.binary, diffF, zF, zNan, binF, args.hdi)
        print(cmd)
        subprocess.run(cmd, shell=True)

        # Generate threshold file
        threshF = "%s/threshFiles/%s_thresh.tsv" % (args.outDir, base)
        io.writeList(["%s\t%s" % (zF, args.zThresh_t1), "%s\t%s" % (args.zscore, args.zThresh_Neg)], threshF, delim="\n")
        
        # Generate rep file lists for input into enrich module
        repF = "%s/nameLists/repFiles/%s.tsv" % (args.outDir, base)
        io.writeList(["\t".join(t1), "\t".join(t2)], repF, delim="\n")
        
        # Run enrich module
        if args.raw:
            cmd = '%s enrich -t %s -s %s -r %s --raw_score_constraint %s -x _%s.txt -f %s_enrichFailReasons.tsv -o %s >> enrich.out' % (args.binary, threshF, repF, args.raw, args.rawThresh, base, base, "%s/enrich" % args.outDir)
        else:
            cmd = '%s enrich -t %s -s %s -x _%s.txt -f %s_enrichFailReasons.tsv -o %s >> enrich.out' % (args.binary, threshF, repF, base, base, "%s/enrich" % args.outDir)

        print(cmd)
        subprocess.run(cmd, shell=True)
        
        # Run deconv
    
#     if matplotReady:
#         if args.raw:
#             #Generate reac counts file
#             args.readCounts = "%s_RC.tsv" % (base)
#             cmd = '%s info -i %s -c %s >> info.out' % (args.binary, args.raw, args.readCounts)
#             print(cmd)
#             subprocess.run(cmd, shell=True)
#         
#             #Read in counts
#             rcD = io.fileDictHeader(args.readCounts, "Sample name", "Sum of probe scores")
#             boxplot(list(rcD.values()), "readCountBoxplot.png", args.rawThresh)
#     
#             if args.thresh and args.pairs and base:
#                 enrFiles = glob.glob("%s/*enriched.txt" % (enrDir))
#                 enrCounts = [len(io.fileList(f, header=False)) for f in enrFiles]
#                 if len(enrCounts) > 0:
#                     boxplot(enrCounts, "enrichedCountBoxplot.png", "200")
# 
#         #Norm read count scatterplot generation
#         if args.repCSscatters and args.colsum and args.pairs:
#             if not os.path.isdir("colsumRepScatters"):
#                 os.mkdir("colsumRepScatters")
#             else:
#                 print("Warning: 'colsumRepScatters' already exists! Col_Sum replicate scatterplots will be placed within existing directory.")
#             
#             pD = io.fileDictLists(args.pairs, header=False)
#             csD = io.fileDictFull(args.colsum, valType="float", rowNames=True)
#             for r1, r2 in pD.items():
#                 for r3 in r2:
#                     scatter(csD[r1], r1, csD[r3], r3, "colsumRepScatters/%s_%s_CS.%s" % (r1, r3, args.scatterFormat), plotLog=True)
# 
#         #Z score scatterplot generation
#         if args.repZscatters and args.zscore and args.pairs:
#             if not os.path.isdir("zRepScatters"):
#                 os.mkdir("zRepScatters")
#             else:
#                 print("Warning: 'zRepScatters' already exists! Z score replicate scatterplots will be placed within existing directory.")
#             
#             pD = io.fileDictLists(args.pairs, header=False)
#             csD = io.fileDictFull(args.zscore, valType="float", rowNames=True)
#             for r1, r2 in pD.items():
#                 for r3 in r2:
#                     scatter(csD[r1], r1, csD[r3], r3, "zRepScatters/%s_%s_Z.%s" % (r1, r3, args.scatterFormat), plotLog=False)
#                 
#         #Enriched peptides scaterplot generation
#         if args.enrichedScatters and args.colsum:
#             cmd = ('%s -d %s -e %s --enrExt _enriched.txt -x NegativeControl -o %s/enrichedScatterplots --plotLog 1' % (args.enrichedScatters, args.colsum, enrDir, enrDir))
#             
#             if args.negNormMat:
#                 cmd += (' --negMatrix %s' % args.negNormMat)
#                 
#             if args.negative_id:
#                 cmd += (' -i %s' % args.negative_id)
#                 
#             if args.negative_names:
#                 cmd += (' -c %s' % args.negative_names)
#                     
#             subprocess.run(cmd, shell=True)
#         elif args.enrichedScatters:
#             print("Warning: plots will not be generated because colsum is not set")


#----------------------End of main()

def sharedSubStr(strL):
    stop=0
    for i in range(min([len(s) for s in strL])):
        if len(set([s[i] for s in strL])) == 1:
            stop+=1
        else:
            return strL[0][:stop]
    return strL[0][:stop]
    
# def scatter(x, xName, y, yName, outname, plotLog=True):
#     #Generate plot
#     fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')            
#     
#     if plotLog:
#         ax.scatter(np.log10(np.array(x)+1), np.log10(np.array(y)+1), alpha=0.5, s=5)
#         ax.set_xlabel("%s log10(score+1)" % xName, fontsize=15)
#         ax.set_ylabel("%s log10(score+1)" % yName, fontsize=15)
#         
#     else:
#         ax.scatter(x, y, alpha=0.5, s=5)
#         ax.set_xlabel("%s" % xName, fontsize=15)
#         ax.set_ylabel("%s" % yName, fontsize=15)
# 
#     fig.savefig(outname, dpi=300, bbox_inches='tight')
#     plt.close(fig)


# def makeDirName(args):
#     dirName = ""
#     infD = io.fileDict(args.thresh, header=False)
#     for k, v in infD.items():
#         typ = k.split("_")[-1].split(".")[-2]
#         val = v.replace(",", "-")
#         dirName+="%s%s_" % (val, typ)
#     if args.raw:
#         dirName+="%sraw" % (args.rawThresh.replace(",", "-"))
#         return dirName
#     else:
#         return dirName[:-1]

# def boxplot(counts, outname, thresh = None):
#     fig,ax = plt.subplots(1,1,figsize=(4, 4),facecolor='w')
#     ax.boxplot([float(x) for x in counts])
#     if thresh:
#         ax.hlines([float(x) for x in thresh.split(",")], 0.5, 1.5, linestyle="--")
#     plt.savefig(outname,dpi=300,bbox_inches='tight')
    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

