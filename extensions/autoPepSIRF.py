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
    inArgs.add_argument("-r", "--raw", help="Input raw count matrix. This is an optional starting place.")
    inArgs.add_argument("-c", "--colsum", help="Input colsum normalized count matrix. This is an optional starting place.")
    inArgs.add_argument("-d", "--diff", help="Input diff normalized count matrix. This is an optional starting place.")
    inArgs.add_argument("-f", "--diffratio", help="Input diff_ratio normalized count matrix. This is an optional starting place.")
    inArgs.add_argument("-b", "--bins", help="Peptide bin file. Required for calculating Z scores.")
    inArgs.add_argument("-z", "--zscore", help="Z score matrix. This is an optional starting place.")
    inArgs.add_argument("-n", "--names", help="File containing sample names. This is an optional starting place.")
    inArgs.add_argument("-p", "--pairs", help="File containing sample pairs. This is an optional starting place.")
    inArgs.add_argument("-t", "--thresh", help="Threshold file for pEnrich. This is an optional starting place.")

    
    threshArgs = parser.add_argument_group('thresholds')
    threshArgs.add_argument("--zThresh", default="6,10", help="Z score threshold. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument("--csThresh", default="20", help="Colum-Sum norm count threshold. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument("--sbdrThresh",  default=False, help="Negative control ratio threshold. Can include up to two floating points separated by a comma.")
    threshArgs.add_argument('--rawThresh',  default="488000", help='Total raw read count for a sample to be included in enrichment analyses. Can include up to two floating points separated by a comma.')
    threshArgs.add_argument('--hdi',  default=0.95, type=float, help='The highest density interval to be used for calculation of mean and stdev in the zscore module.')

    controlArgs = parser.add_argument_group('control info')
    controlArgs.add_argument("--negNormMat", help="Alternative colsum normalized matrix form which to obtain the negative controls.")
    controlArgs.add_argument("--negative_id", help="Optional approach for identifying negative controls. Provide a unique string at the start of all negative control samples.")
    controlArgs.add_argument("--negative_names", help="Optional approach for identifying negative controls. Comma-separated list of negative control sample names.")

    enrichArgs = parser.add_argument_group('enrich options')
    enrichArgs.add_argument("--sEnrich", default=False, action="store_true", help="Generate lists of enriched peptides separately for each pulldown. Will actually run p_enrich, but with the same sample specified for each replicate.")
    enrichArgs.add_argument("--inferPairs", default=False, action="store_true", help="Infer sample pairs from names. This option assumes names of replicates will be identical with the exception of a final string denoted with a '_'. For example, these names would be considered two replicates of the same sample: VW_100_1X_A and VW_100_1X_B")

    outArgs = parser.add_argument_group('output options')
    enrichArgs.add_argument("--repZscatters", default=False, action="store_true", help="Generate scatter plots comparing Z scores for all sample pairs. Even if set at command line, will be turned off when '--sEnrich' is used.")
    enrichArgs.add_argument("--repCSscatters", default=False, action="store_true", help="Generate scatter plots comparing colu_sum normalized scores for all sample pairs. Even if set at command line, will be turned off when '--sEnrich' is used.")
    enrichArgs.add_argument("--scatterFormat", default="png", help="Output file format for replicate scatterplots.")
    enrichArgs.add_argument("--enrichedScatters", default=False, help="Generate scatter plots comparing col-sum normalized read counts between a sample and negative controls. Argument provided should be the path to enrichedScatterplots.py")

    args = parser.parse_args()
    
    #Create base string for output files
    if args.raw:
        base = ".".join(args.raw.split(".")[:-1])
    elif args.colsum:
        base = ".".join(args.colsum.split(".")[:-1])[:-3]
    elif args.diff:
        base = ".".join(args.diff.split(".")[:-1])[:-4]
    elif args.diffratio:
        base = ".".join(args.diffratio.split(".")[:-1])[:-5]
    elif args.zscore:
        base = ".".join(args.zscore.split(".")[:-1])[:-2]
    else:
        base = None
        print("Not going to be able to do much without a score matrix of some type.")
    
    
    # If a raw count matrix is provided, but a colusum norm matrix is NOT
    if args.raw and not args.colsum:
        args.colsum = "%s_CS.tsv" % (base)
        cmd = "%s norm -a col_sum -p %s -o %s >> norm.out" % (args.binary, args.raw, args.colsum)
        print(cmd)
        subprocess.run(cmd, shell=True)
#        subprocess.run([args.binary, "norm -a col_sum -p", args.raw, "-o", args.colsum, ">> norm.out"])

    if args.colsum:
        
        if not args.negative_id and not args.negative_names:
            print("You must proivde either '--negative_id' or '--negative_names' to allow negative control based normalization.")

        else:

            # Generate string for specifying how to get neg control data
            negInfo = ""
            if args.negNormMat:
                negInfo+=' --negative_control "%s" ' % (args.negNormMat)
            if args.negative_id:
                negInfo+=" --negative_id %s " % (args.negative_id)
            if args.negative_names:
                negInfo+=" --negative_names %s " % (args.negative_names)

            # Generate other normalized files
            if not args.diff:
                args.diff = "%s_SBD.tsv" % (base)
               
                cmd = "%s norm -a diff -p %s -o %s %s >> norm.out" % (args.binary, args.colsum, args.diff, negInfo)
                print(cmd)
                subprocess.run(cmd, shell=True)
#                subprocess.run([opts.binary, "norm -a diff -p", args.colsum, "-o", args.diff, negInfo, ">> norm.out"])

            if not args.diffratio:
                args.diffratio = "%s_SBDR.tsv" % (base)

                cmd = "%s norm -a diff_ratio -p %s -o %s %s >> norm.out" % (args.binary, args.colsum, args.diffratio, negInfo)
                print(cmd)
                subprocess.run(cmd, shell=True)
#                subprocess.run([opts.binary, "norm -a diff_ratio -p", args.colsum, "-o", args.diffratio, negInfo, ">> norm.out"])
        
    if args.bins and args.diff:
        # Generate Z scores
        args.zscore = "%s_Z-HDI%d.tsv" % (base, int(args.hdi*100))
        args.znan = "%s_Z-HDI%d.nan" % (base, int(args.hdi*100))
        cmd = '%s zscore -s %s -o %s -n %s -b "%s" -d %f >> zscore.out' % (args.binary, args.diff, args.zscore, args.znan, args.bins, args.hdi)
        print(cmd)
        subprocess.run(cmd, shell=True)
        
    elif not args.bins:
        print("You must proivde '--bins' for Z score calculation.")

    # Generate list of sample names, if not provided
    if not args.names and base:
        args.names = "%s_SN.tsv" % (base)
        if args.raw:
            cmd = '%s info -i %s -s %s >> info.out' % (args.binary, args.raw, args.names)
            print(cmd)
            subprocess.run(cmd, shell=True)
        elif args.colsum:
            cmd = '%s info -i %s -s %s >> info.out' % (args.binary, args.colsum, args.names)
            print(cmd)
            subprocess.run(cmd, shell=True)
        elif args.zscore:
            cmd = '%s info -i %s -s %s >> info.out' % (args.binary, args.zscore, args.names)
            print(cmd)
            subprocess.run(cmd, shell=True)
        elif args.diff:
            cmd = '%s info -i %s -s %s >> info.out' % (args.binary, args.diff, args.names)
            print(cmd)
            subprocess.run(cmd, shell=True)
        elif args.diffratio:
            cmd = '%s info -i %s -s %s >> info.out' % (args.binary, args.diffratio, args.names)
            print(cmd)
            subprocess.run(cmd, shell=True)
        else:
            print("No file was providing for generating a list of sample names.")
    
    
    # Generate pairs file, if not provided
    if not args.pairs and args.names and base:
        sNames = io.fileList(args.names, header=False)

        if args.sEnrich and args.names:
            args.pairs = "%s_pseudoPN.tsv" % (base)
            with open(args.pairs, "w") as fout:
                for sn in sNames:
                    fout.write("%s\t%s\n" % (sn, sn))
            
            # Turn off replicate scatterplot generation, as these will not be true replicates
            args.repZscatters = False
            args.repCSscatters = False
            
        elif args.inferPairs:
            pDict = defaultdict(list)
            for each in sNames:
                simple = "_".join(each.split("_")[:-1])
                pDict[simple].append(each)
                
            args.pairs = "%s_PN.tsv" % (base)
            with open(args.pairs, "w") as fout:
                for k, v in pDict.items():
                    if len(v) == 2:                       # If exactly two replicates were found, add them to the pairs file
                        fout.write("%s\n" % ("\t".join(v)))
                    elif len(v) > 2:                       # If more than two replicates were found, add all possible pairs
                        for a,b in it.combinations(v,2):
                            fout.write("%s\t%s\n" % (a, b))
                    else:
                        print("Only one replicate found for %s: %s" % (simple, each))

        else:
            print("To run p_enrich module, you must provide one of the following: '--pairs', '--inferPairs', '--sEnrich'")

    # Generate threshold file
    if not args.thresh and base:
        args.thresh = "%s_thresh.tsv" % (base)
        with open(args.thresh, "w") as fout:
            if args.zscore:
                fout.write("%s\t%s\n" % (args.zscore, args.zThresh))
            if args.colsum:
                fout.write("%s\t%s\n" % (args.colsum, args.csThresh))
            if args.diffratio and args.sbdrThresh:
                fout.write("%s\t%s\n" % (args.diffratio, args.sbdrThresh))
    
    # Run p_enrich module
    if args.thresh and args.pairs and base:
        enrDir = makeDirName(args)
        if args.raw:
            cmd = '%s p_enrich -t %s -s %s -r %s --raw_score_constraint %s -x _enriched.txt -o %s >> penrich.out' % (args.binary, args.thresh, args.pairs, args.raw, args.rawThresh, enrDir)
        else:
            cmd = '%s p_enrich -t %s -s %s -x _enriched.txt -o %s >> penrich.out' % (args.binary, args.thresh, args.pairs, enrDir)

        print(cmd)
        subprocess.run(cmd, shell=True)
    # Turn off enriched scatterplots generation, as they cannot be run without a list of enriched peptides
    elif args.enrichedScatterplots:
        print("Warning: p_enrich module not run and list of enriched peptides not generated. --enrichedScatterplots will be set to False.")
        args.enrichedScatterplots = False

    
    if matplotReady:
        if args.raw:
            #Generate reac counts file
            args.readCounts = "%s_RC.tsv" % (base)
            cmd = '%s info -i %s -c %s >> info.out' % (args.binary, args.raw, args.readCounts)
            print(cmd)
            subprocess.run(cmd, shell=True)
        
            #Read in counts
            rcD = io.fileDictHeader(args.readCounts, "Sample name", "Sum of probe scores")
            boxplot(list(rcD.values()), "readCountBoxplot.png", args.rawThresh)
    
            if args.thresh and args.pairs and base:
                enrFiles = glob.glob("%s/*enriched.txt" % (enrDir))
                enrCounts = [len(io.fileList(f, header=False)) for f in enrFiles]
                if len(enrCounts) > 0:
                    boxplot(enrCounts, "enrichedCountBoxplot.png", "200")

        #Norm read count scatterplot generation
        if args.repCSscatters and args.colsum and args.pairs:
            if not os.path.isdir("colsumRepScatters"):
                os.mkdir("colsumRepScatters")
            else:
                print("Warning: 'colsumRepScatters' already exists! Col_Sum replicate scatterplots will be placed within existing directory.")
            
            pD = io.fileDict(args.pairs, header=False)
            csD = io.fileDictFull(args.colsum, valType="float", rowNames=True)
            for r1, r2 in pD.items():
                scatter(csD[r1], r1, csD[r2], r2, "colsumRepScatters/%s_%s_CS.%s" % (r1, r2, args.scatterFormat), plotLog=True)

        #Z score scatterplot generation
        if args.repZscatters and args.zscore and args.pairs:
            if not os.path.isdir("zRepScatters"):
                os.mkdir("zRepScatters")
            else:
                print("Warning: 'zRepScatters' already exists! Z score replicate scatterplots will be placed within existing directory.")
            
            pD = io.fileDict(args.pairs, header=False)
            csD = io.fileDictFull(args.zscore, valType="float", rowNames=True)
            for r1, r2 in pD.items():
                scatter(csD[r1], r1, csD[r2], r2, "zRepScatters/%s_%s_Z.%s" % (r1, r2, args.scatterFormat), plotLog=False)
                
        #Enriched peptides scaterplot generation
        if args.enrichedScatters and args.colsum:
            cmd = ('python3 %s -d %s -e %s --enrExt _enriched.txt -x NegativeControl -o %s/enrichedScatterplots --plotLog 1' % (args.enrichedScatters, args.colsum, enrDir, enrDir))
            
            if args.negNormMat:
                cmd += (' --negMatrix %s' % args.negNormMat)
                
            if args.negative_id:
                cmd += (' -i %s' % args.negative_id)
                
            if args.negative_names:
                cmd += (' -c %s' % args.negative_names)
                    
            subprocess.run(cmd, shell=True)
        elif args.enrichedScatters:
            print("Warning: plots will not be generated because colsum is not set")


#----------------------End of main()

def scatter(x, xName, y, yName, outname, plotLog=True):
    #Generate plot
    fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')            
    
    if plotLog:
        ax.scatter(np.log10(np.array(x)+1), np.log10(np.array(y)+1), alpha=0.5, s=5)
        ax.set_xlabel("%s log10(score+1)" % xName, fontsize=15)
        ax.set_ylabel("%s log10(score+1)" % yName, fontsize=15)
        
    else:
        ax.scatter(x, y, alpha=0.5, s=5)
        ax.set_xlabel("%s" % xName, fontsize=15)
        ax.set_ylabel("%s" % yName, fontsize=15)

    fig.savefig(outname, dpi=300, bbox_inches='tight')
    plt.close(fig)


def makeDirName(args):
    dirName = ""
    infD = io.fileDict(args.thresh, header=False)
    for k, v in infD.items():
        typ = k.split("_")[-1].split(".")[-2]
        val = v.replace(",", "-")
        dirName+="%s%s_" % (val, typ)
    if args.raw:
        dirName+="%sraw" % (args.rawThresh.replace(",", "-"))
        return dirName
    else:
        return dirName[:-1]

def boxplot(counts, outname, thresh = None):
    fig,ax = plt.subplots(1,1,figsize=(4, 4),facecolor='w')
    ax.boxplot([float(x) for x in counts])
    if thresh:
        ax.hlines([float(x) for x in thresh.split(",")], 0.5, 1.5, linestyle="--")
    plt.savefig(outname,dpi=300,bbox_inches='tight')
    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

