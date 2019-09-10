#!/usr/bin/env python

import optparse, re

# Used to check the taxonomy categories assigned to protein sequences
#      Based on kmer overlap to validated reference sequences

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-r', '--ref', help='Fasta file with validated reference protein sequences. [None, REQ]')
    p.add_option('-q', '--queries', help='Fasta file with query protein sequences. [None, REQ]')
    p.add_option('-o', '--out', help='Base string for output files. [None, REQ]')
    p.add_option('-k', '--kmer', type='int', default=4, help='Kmer size. [4]')
    p.add_option('-t', '--thresh', type='float', default=0.2, help='Minimum proportion of shared kmers for a protein to be considered a potential member of one of the refernce taxa. [0.2]')
    p.add_option('-c', '--cat', default="sp", help='Taxon category of interest. Must be one of these: "str", "sp", "gen", "fam" [sp]')
    opts, args = p.parse_args()
    
    refD = read_fasta_dict_upper(opts.ref)        #Read in reference seqs
    refTaxKmersD = kmersByTaxon(refD, opts)       #Arrange by taxon level of interest and generate kmers
    
    queD = read_fasta_dict_upper(opts.queries)              #Read in query seqs
    putAssign = checkQueries(queD, refTaxKmersD, opts)      #Compare queries to references
    checkAssign(putAssign, queD, opts)
#----------------------End of main()

def checkAssign(putAssign, queD, opts):
    goodNames=[]
    goodSeqs=[]
    badNames=[]
    badSeqs=[]
    
    with open("%s_good.txt" % opts.out, "w") as foutGood, open("%s_bad.txt" % opts.out, "w") as foutBad:
        foutGood.write("ProtName\tProtLength\tID\tBestMatch\tProp%dmers\n" % opts.kmer)
#        foutGood.write("ProtName\tProtLength\tMatchID\tBestMatch\tBestProp%dmers\tProtID\tBestProtIDMatch\tProp%dmers" % opts.kmer)
        foutBad.write("ProtName\tProtLength\tProtID\tMatchID\tBestMatch\tProp%dmers\n" % opts.kmer)
        for matchID, info in putAssign.items():
            for name, other in info.items():
                prop, rN = other
                thisID = parseTax(name, opts.cat)
                if thisID == matchID:
                    foutGood.write("%s\t%d\t%s\t%s\t%.3f\n" % (name, len(queD[name]), thisID, rN, prop))
                    goodNames.append(name)
                    goodSeqs.append(queD[name])
                else:
                    foutBad.write("%s\t%d\t%s\t%s\t%s\t%.3f\n" % (name, len(queD[name]), thisID, matchID, rN, prop))
                    badNames.append(name)
                    badSeqs.append(queD[name])
    if goodNames:
        write_fasta(goodNames, goodSeqs, "%s_good.fasta" % (opts.out))
    if badNames:
        write_fasta(badNames, badSeqs, "%s_bad.fasta" % (opts.out))


def checkQueries(queD, refTaxKmersD, opts):
    outNames=[]
    outSeqs=[]
    with open("%s_missed.txt" % opts.out, "w") as foutMissed:
        foutMissed.write("ProtName\tProtLength\tProtID\tBestMatch\tProp%dmers\n" % opts.kmer)
        putAssign = {x:{} for x in refTaxKmersD}
        for name, seq in queD.items():
            bestHits = {x:[0,""] for x in refTaxKmersD}
            ks = kmers(seq, opts.kmer)
            for cat, info in refTaxKmersD.items():
                for rN, rK in info.items():
                    ovlp = compPairKs(ks,rK)
                    if ovlp>bestHits[cat][0]:
                        bestHits[cat] = [ovlp,rN]
            maxVal, maxName, maxTax = findMax(bestHits)
            if maxVal>=opts.thresh:
                putAssign[maxTax][name] = [maxVal, maxName]
            else:                                            #Check to see if seqs without good matches are assigned to focal taxa
                thisID = parseTax(name, opts.cat)
                if thisID in refTaxKmersD:
                    foutMissed.write("%s\t%d\t%s\t%s\t%.3f\n" % (name, len(queD[name]), thisID, maxName, maxVal))
                    outNames.append(name)
                    outSeqs.append(queD[name])
    if outNames:
        write_fasta(outNames, outSeqs, "%s_missed.fasta" % (opts.out))
    return putAssign

def findMax(hitD):
    mx = 0
    info=[0, "", ""]
    for t, i in hitD.items():
        if i[0]>mx:
            mx=i[0]
            info=[i[0],i[1],t]
    return info

def compPairKs(aK,bK):
    ovlp = aK.intersection(bK)
    if len(aK)<=len(bK):
        return len(ovlp)/len(aK)
    else:
        return len(ovlp)/len(bK)


def kmersByTaxon(refD, opts):
    kD = {}
    for name, seq in refD.items():
        id = parseTax(name, opts.cat)
        if id not in kD:
            kD[id]={}
        kD[id][name] = kmers(seq, opts.kmer)
    return kD

def kmers(seq,k, step=1):
    out=[]
    for i in range(0, len(seq)-k+1, step):
        out.append(seq[i:i+k])
    return set(out)

def parseTax(name, category):
    catD = {"str":1,"sp":2,"gen":3,"fam":4}
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_id = oxpat.search(name)
    if tax_id:
        return tax_id.group(catD[category])
    else:
        #print(name)
        return None

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

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

