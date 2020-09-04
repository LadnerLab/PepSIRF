#!/usr/bin/env python

import optparse
from scipy import stats

#This script takes as input a file containing normalized count data
    #As well as a file defining group assignments
    #And performs T-tests for each probe

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-d', '--data', help='Count data file. [None, REQ]')
    p.add_option('-g', '--groups', help='Tab delimited file with two rows, containing the names of samples in each group. One group per row. [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file [None, REQ]')
    p.add_option('--minMax', default=100, type='int', help='At least one sample, in one group needs to have a count >= this value for a test to be run [100]')    
    p.add_option('--equalVar', default=False, action="store_true", help='Use this flag if you want to assume equal population variances. [False]')

    
    opts, args = p.parse_args()
    
    #Read in groups
    grp1, grp2 = readGroups(opts.groups)
    print(grp1, grp2)
    #Read in data
    data = parseCounts(opts.data)
    
    with open(opts.out, "w") as fout:
        fout.write("Probe\tT-statistic\tTwo tailed P-value\n")
        #Step through probes and calculate p-values
        tests=0
        for p in data[grp1[0]].keys():
            one = [data[x][p] for x in grp1]
            two = [data[x][p] for x in grp2]
            if max(one+two) >= opts.minMax:
                tests+=1
                T, pval = stats.ttest_ind(one,two, equal_var = opts.equalVar)
                fout.write("%s\t%.3f\t%.9f\n" % (p,T,pval))
    print("%d Tests were run, in total" % (tests))
#----------------------End of main()

def readGroups(grpFile):
    groups = []
    with open(grpFile, "r") as fin:
        for line in fin:
            cols = line.rstrip().split("\t")
            if cols:
                groups.append(cols)
    return groups[0], groups[1]

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


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

