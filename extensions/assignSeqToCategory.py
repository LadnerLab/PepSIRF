#!/usr/bin/env python

import optparse, re
import inout as io
import fastatools as ft
import kmertools as kt
from collections import defaultdict

# Adapted from checkTaxonomy.py to be a more general kmer based approach for assigning 
#      sequenced to categories through comparison with validated reference sequences

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-r', '--ref', help='Fasta file with validated reference protein sequences. [None, REQ]')
    p.add_option('-q', '--queries', help='Fasta file with query protein sequences. [None, REQ]')
    p.add_option('-o', '--out', help='Base string for output files. [None, REQ]')
    p.add_option('-k', '--kmer', type='int', default=4, help='Kmer size. [4]')
    p.add_option('-t', '--thresh', type='float', default=0.2, help='Minimum proportion of shared kmers for a protein to be considered a potential member of one of the refernce taxa. [0.2]')
    p.add_option('-m', '--meta', help='Tab-delim metadata file with headers. One column should correspond to simple names of the reference fastas. Another should correspond to the category of interest.  [None, REQ]')
    p.add_option('-c', '--cat', help='Metadata header for category of interest [None, REQ]')
    p.add_option('-n', '--names', help='Metadata header for seq names. It is assumed that these are "simple" names, corresponding to the fasta file names prior to the first whitespace character [None, REQ]')
    p.add_option('--minRefLen', type="int", default=40, help='Minimum required length (in AA) for a reference sequence to be used [40]')

    opts, args = p.parse_args()
    
    # Read in reference sequences
    refD = ft.read_fasta_dict_upper(opts.ref)        #Read in reference seqs
    simpleRefD = {n.split()[0]:s for n,s in refD.items() if len(s)>=opts.minRefLen}
    
    # Read in data of interest from metadata file
    catD = io.fileDictHeader(opts.meta, opts.names, opts.cat)
    
    #Read in query sequences
    queD = ft.read_fasta_dict_upper(opts.queries)              #Read in query seqs

    #Compare queries to references, generate output files
    assign(simpleRefD, queD, catD, opts)      
#----------------------End of main()



def assign(refD, queD, catD, opts):
    allCats = sorted(list(set([catD[n] for n in refD])))

    outD={c:[[],[]] for c in allCats}

    with open("%s_hitInfo.tsv" % opts.out, "w") as fout:
        fout.write("SeqName\tSeqLength\tBestCat\t%s\n" % ("\t".join(allCats)))

        for qName, qSeq in queD.items():
            theseHits = {}
        
            for rName, rSeq in refD.items():
                cat = catD[rName]
                ovlp = kt.compSeqs(qSeq, rSeq, opts.kmer)
                if cat not in theseHits:
                    theseHits[cat] = (ovlp, len(rSeq))
                elif ovlp>theseHits[cat][0]:
                    theseHits[cat] = (ovlp, len(rSeq))
        
            bestCat, maxProp = findMax(theseHits)
        
            fout.write("%s\t%d\t%s\t%s\n" % (qName, len(qSeq), bestCat, "\t".join(["%.3f" % (theseHits[c][0]) for c in allCats])))

            if maxProp>=opts.thresh:
                outD[bestCat][0].append(qName)
                outD[bestCat][1].append(qSeq)
    
    for cat, infL in outD.items():
        if infL[0]:
            ft.write_fasta(infL[0], infL[1], "%s_%s.fasta" % (opts.out, cat))


def findMax(hitD):
    mx = 0
    mxCat=""
    mxRefLen=0
    for c, infoT in hitD.items():
        p,l = infoT
        if p>mx:
            mx=p
            mxCat=c
            mxRefLen=l
        elif p==mx and l<mxRefLen:
            mx=p
            mxCat=c
            mxRefLen=l

    return mxCat, mx



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

