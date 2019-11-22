#!/usr/bin/env python

import optparse

#Code or decode probe names in various file types

def main():
    usage = '%prog [options] probes1.txt [probes2.txt ...]'
    p = optparse.OptionParser()

    p.add_option('-m', '--map',  help='File containing map to link names in enriched list to those in pep fasta. [None, REQ]')
    p.add_option('-f', '--fastas',  help='One or more fasta files, comma separated. Expect that the seq names are probe names, either coded or not. [OPT]')
    p.add_option('-l', '--links',  help='One or more linkage maps, comma separated. [OPT]')
    p.add_option('-s', '--lists',  help='One or more probe lists, comma separated. [OPT]')

    opts, args = p.parse_args()

    # Read in code map
    codeD = fileDict(opts.map)
    revCodeD = {v:k for k,v in codeD.items()}

    if opts.links:
        for each in opts.links.split(","):
            eachD = fileDict(each, header=True)
            with open("%s_coded.tsv" % each[:-4], "w") as fout:
                fout.write("Name\tSpecies\n")
                for p,spec in eachD.items():
                    if p in codeD:
                        fout.write("%s\t%s\n" % (codeD[p], spec))
                    elif p in revCodeD:
                        fout.write("%s\t%s\n" % (revCodeD[p], spec))
                    else:
                        "Don't recognize this name: %s" % p

#----------------------End of main()

def fileDict(file, header=False):
    D={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if not header or lc>1:
                cols = line.rstrip("\n").split("\t")
                D[cols[0]] = cols[1]
    return D


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
