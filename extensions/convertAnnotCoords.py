#!/usr/bin/env python

import optparse, re

#Used to take annotations for an individual sequence and translate them to an alignment including that sequence

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--info', help='File with info on sequence, coordinates, alignments. Coordinates should be 1-based. [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file. [None, REQ]')
    opts, args = p.parse_args()
    
    with open(opts.out, "w") as fout:
        #Step through each line in the info file
        with open(opts.info, "r") as fin:
            lc=0
            for line in fin:
                lc+=1
                if lc>1:
                    cols = line.strip("\n").split("\t")
                    fD = read_fasta_dict_upper(cols[4])
                    seq2alD = mapAlPos(fD[cols[1]])
                    newAnnots = convertAnnots(cols[2], seq2alD)
                    fout.write("%s\t%s\t%s\t%s\t%s\n" % (cols[0], cols[1], newAnnots, cols[3], cols[4]))
                else:
                    fout.write(line)

#----------------------End of main()

def convertAnnots(start, map):
    info = [x.split(",") for x in start.split("~")]
    new = ["%s,%s,%s" % (map[x[0]],map[x[1]],x[2]) for x in info]
    return "~".join(new)

def mapAlPos(seq):
    map={}
    a=0
    s=0
    for pos in seq:
        a+=1
        if pos != "-":
            s+=1
            map[str(s)]=str(a)
    return map

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



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

