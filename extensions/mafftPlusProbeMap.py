#!/usr/bin/env python

import optparse, subprocess, re

#This script takes as input a set of fasta files, aligns each using Mafft and then generates an info file for probe-level alignments

def main():
    usage = '%prog [options] fasta1  [fasta2 ...]'
    p = optparse.OptionParser()
    p.add_option('-p', '--probes', default="/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/peptides/unrepSW_9_30_70_design_combined_wControls.fasta_unique", help='Fasta with probe seqs. [/Volumes/GoogleDrive/Shared drives/MyImmunity/PanviralDesign/PV1/peptides/unrepSW_9_30_70_design_combined_wControls.fasta_unique]')
    p.add_option('--mafft', default="mafft", help='How to call Mafft aligner on your machine. [mafft]')
    p.add_option('-i', '--iter', default=5,  type='int', help='Max # of iterations for Mafft [5]')
    p.add_option('--skipMafft', default=False,  action='store_true', help='Use this flag if the provided fasta files are already aligned [False]')

    opts, args = p.parse_args()
    
    for fasta in args:
        if not opts.skipMafft:
            alFasta = "%s_mafft-maxiter%d.fasta" % (".".join(fasta.split(".")[:-1]), opts.iter)
            subprocess.run([opts.mafft, "--quiet", "--maxiterate", opts.iter, fasta, ">", alFasta])
        else:
            alFasta = fasta
        
        #Make dict that links target seq position to alignment positions (1-based)
        alignD = parseAlignment(alFasta)
    
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
                    nonIndel = []
                    for p in range(start, stop+1):
                        if p in alignD[tarName]:
                            nonIndel.append(str(p))
                    fout.write("%s\t%d\t%d\t%s\n" % (each,start,stop, "~".join(nonIndel)))

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

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

