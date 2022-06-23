#!/usr/bin/env python

import optparse

#This script takes as input a set of aligned fasta files and generates an info file for probe-level alignments

def main():
    usage = '%prog [options] fasta1  [fasta2 ...]'
    p = optparse.OptionParser()
    p.add_option('-p', '--probes', default="/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/library_design/oligos/2019-02-27/unrepSW_9_30_70_design_combined_wControls_2019-09-12.fasta_unique", help='Fasta with probe seqs. [/Volumes/GoogleDrive/Shared drives/LadnerLab/Projects/panviral_pepseq/library_design/oligos/2019-02-27/unrepSW_9_30_70_design_combined_wControls_2019-09-12.fasta_unique]')
    p.add_option('--subsetAlignments', default=False, action="store_true", help='Use this option if the seqs in the alignment are NOT full length target sequences. In this case, the sequence names should end with a string such as " 762-962" indicating the positions of the extracted region in X-coordinates. [False]')

    opts, args = p.parse_args()
    
    for alFasta in args:
        #Make dict that links target seq position to alignment positions (1-based)
        alignD = parseAlignment(alFasta)
        
        #If seqs in alignment are subsets of the original target sequences, then alter alignD and generate a supporting data object
        if opts.subsetAlignments:
            name2coords = {" ".join(n.split(" ")[:-1]): [int(x) for x in n.split(" ")[-1].split("-")] for n in alignD}
            alignD = {" ".join(n.split(" ")[:-1]): s for n,s in alignD.items()}
            for k, v in alignD.items():
                alignD[k] = {p+name2coords[k][0]-1:a for p,a in v.items()}
        
        #Extract appropriate probes locations and write to out file
        outFile = "%s_probesAligned.txt" % (".".join(alFasta.split(".")[:-1]))
        with open(outFile, "w") as fout:
            fout.write("ProbeName\tAlignStart(1-based)\tAlignStop(1-based)\tAlignPos\n")
            #Read in probe seqs
            pNames, pSeqs = read_fasta_lists(opts.probes)
            #Find probes designed from seqs in alignment and extract alignment locations
            for i, each in enumerate(pNames):
                tarName, start, stop = parseProbeName(each)
                if tarName in alignD:
                    if start in alignD[tarName] and stop in alignD[tarName]:
                        nonIndel = []
                        for p in range(start, stop+1):
                            if p in alignD[tarName]:
                                nonIndel.append(str(alignD[tarName][p]))
    #                            nonIndel.append(str(p))
                        fout.write("%s\t%d\t%d\t%s\n" % (each, alignD[tarName][start], alignD[tarName][stop], "~".join(nonIndel)))
#                    fout.write("%s\t%d\t%d\t%s\n" % (each,start,stop, "~".join(nonIndel)))

#----------------------End of main()

def parseProbeName(each):
    if each.endswith("fasta"):
        each = each.split()[0]
    parts = each.split("_")
    return "_".join(parts[:-2]), int(parts[-2])+1, int(parts[-1])

def parseAlignment(fsta):
    names, seqs = read_fasta_lists(fsta)
    alD={n:{} for n in names}

    for i,each in enumerate(seqs):
        pos = 0
        for j,aa in enumerate(each):
            if aa != "-":
                pos+=1
                alD[names[i]][pos] = j+1
    return(alD)

    
# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
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

