#!/usr/bin/env python

import optparse, os, subprocess, re

#This script reads in a score matrix and outputs sets of "enriched" probes designed from a given taxon

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-s', '--matrix',  help='Score matrix. [None, REQ]')
    p.add_option('-m', '--map',  help='File containing map to link coded names to names with tax info. [None, REQ]')
    p.add_option('--thresh', default=10, type='int', help='Inclusive threshold for output to file [10]')
    p.add_option('--mapOrder', default=1,  type='int', help='Integer indicate the column # of the name labels in the matrix. Must be 1 or 2 (1-based)[1]')
    p.add_option('-i', '--id',  help='Species ID of interest. [REQ]')
    p.add_option('--comboOnly', default=False, action="store_true", help='If used, will only output a single list of probes enriched across all samples. [False]')

    opts, args = p.parse_args()
    
    
    #Read in Score Matrix
    scoreD = parseCounts(opts.matrix)

    #Read in name map, if provided
    if opts.map:
        pepMap = readmap(opts.map, order=opts.mapOrder-1)

    total=[]
    for samp, inf in scoreD.items():
        enr=[]
        for p,s in inf.items():
            if s>=opts.thresh:
                taxInf = parse_tax(pepMap[p])
                if taxInf:
                    if taxInf[2] == opts.id:
                        enr.append(p)
        if not opts.comboOnly:
            with open("%s_id%s_thresh%d.txt" % (samp,opts.id, opts.thresh), "w") as fout:
                fout.write("%s\n" % ("\n".join(enr)))
        total+=enr
    total=list(set(total))
    with open("total_id%s_thresh%d.txt" % (opts.id, opts.thresh), "w") as fout:
        fout.write("%s\n" % ("\n".join(total)))

#----------------------End of main()

def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

def filelist(file):
    l=[]
    with open(file, "r") as fin:
        for line in fin:
            l.append(line.strip("\n"))
    return l

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


#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def parse_tax(name):
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_id = oxpat.search(name)
    if tax_id:
        species,genus,family = (tax_id.group(2),tax_id.group(3),tax_id.group(4))
        return family,genus,species
    else:
        #print(name)
        return None

def speciesNames(file):
    id2name={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1:
                cols=line.strip("\n").split("\t")
                id2name[cols[4]]=cols[0]
    return(id2name)

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
