#!/usr/bin/env python

import optparse, os, subprocess, re

#This script clusters and aligned enriched peptides that share epitopes

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-e', '--enriched',  help='List of enriched peptide names, one per line. [None, REQ]')
    p.add_option('-p', '--peps',  help='Fasta file peptide sequences. [None, REQ]')
    p.add_option('-o', "--outDir", default="epiClusts", help='Name for opts.outDirctory in which output files will be generated. This opts.outDirctory should NOT exist already. [epiClusts]')
    p.add_option('-m', '--map',  help='File containing map to link names in enriched list to those in pep fasta. [OPT]')
    p.add_option('--mapOrder', default=1,  type='int', help='Integer indicate the column # of the name labels in the enriched list. Must be 1 or 2 (1-based)[1]')
    p.add_option('-s', '--species',  help='File containing link between species taxonomy IDs and full names. [OPT]')
    p.add_option('-k', '--ksize', default=6, type='int', help='Kmer size to use to cluster peptides. [6]')
    p.add_option('--minK', default=4, type='int', help='Minimum # of shared kmers to cluster [4]')
    p.add_option('--minProp', default=0.2, type='float', help='Minimum proportion of pairwise comparisons to meet threshold in order to combine clusters. Must be between 0-1. [0.2]')
    p.add_option('--muscle', default="muscle3.8.31_i86darwin64", help='How to call Muscle aligner on your machine. [muscle3.8.31_i86darwin64]')

    opts, args = p.parse_args()
    
    #Read in peptides
    pepDict = read_fasta_dict_upper(opts.peps)

    #Read in Enriched peptides
    peps = filelist(opts.enriched)

    #Read in name map, if provided
    if opts.map:
        pepMap = readmap(opts.map, order=opts.mapOrder-1)
        
    #Read in species ID map, if provided
    if opts.species:
        id2name = speciesNames(opts.species)
        
    #Generate clusters
    if opts.map:
        clu = makeClusts({k:pepDict[pepMap[k]] for k in peps}, opts.ksize, opts.minK)
    else:
        clu = makeClusts({k:pepDict[k] for k in peps}, opts.ksize, opts.minK)
    
    #Merge similar clusters
    cMerged, singles = mergeClusts(clu, opts.ksize, opts.minK, opts.minProp)
    
    #Write out unaligned fasta files and cluster stats
    os.mkdir(opts.outDir)

    with open("%s/clusterInfo.txt" % opts.outDir, "w") as fout:
        #Write out headers for info file
        if opts.map:
            fout.write("Cluster\tNumPeps\tSpIDs\tSpNames\tPepIDs\tPepIDsCoded\n")
        else:
            fout.write("Cluster\tNumPeps\tSpIDs\tSpNames\tPepIDs\n")
            
        outnames=[]
        for i, each in enumerate(cMerged):
            names=[]
            seqs=[]
            spIDs=set()
            for k,v in each.items():
                names.append(k)
                seqs.append(v)
                #Parse taxonomy
                if opts.map:
                    taxInfo = parse_tax(pepMap[k])
                else:
                    taxInfo = parse_tax(k)
                if taxInfo:
                    spIDs.add(taxInfo[2])

            oname = "%s/epiclust_%03d.fasta" % (opts.outDir,i)
            outnames.append(oname)
            write_fasta(names, seqs, oname)
            
            
            #Write out cluster info
            spIDs = sorted(list(spIDs))
            if opts.species:
                spNames = [id2name[x] for x in spIDs]
            else:
                spNames = []
            if opts.map:
                fout.write("%03d\t%d\t%s\t%s\t%s\t%s\n" % (i, len(each), "~".join(spIDs), "~".join(spNames), "~".join([pepMap[k] for k in each.keys()]), "~".join(list(each.keys()))))
            else:
                fout.write("%03d\t%d\t%s\t%s\t%s\n" % (i, len(each), "~".join(spIDs), "~".join(spNames), "~".join(list(each.keys()))))

    #Write out singles together
    with open("%s/singlesInfo.txt" % opts.outDir, "w") as fout:
        if opts.map:
            fout.write("PepIDCoded\tSpID\tSpName\tPepID\n")
        else:
            fout.write("PepID\tSpID\tSpName\n")
        names=[]
        seqs=[]
        for each in singles:
            k=list(each.keys())[0]
            v=list(each.values())[0]
            names.append(k)
            seqs.append(v)
            #Parse taxonomy
            if opts.map:
                taxInfo = parse_tax(pepMap[k])
            else:
                taxInfo = parse_tax(k)
            if taxInfo:
                spID = (taxInfo[2])
            else:
                spID=""
            if opts.species:
                spName = id2name[spID]
            else:
                spName = ""
                
            if opts.map:
                fout.write("%s\t%s\t%s\t%s\n" % (k, spID, spName, pepMap[k]))
            else:
                fout.write("%s\t%s\t%s\n" % (k, spID, spName))

        write_fasta(names, seqs, "%s/unclustered.fasta" % (opts.outDir))
    
    #Generate aligned fasta files 

    alDir = "%s/aligned" % opts.outDir
    os.mkdir(alDir)

    for each in outnames:
        subprocess.run([opts.muscle, "-in", each, "-out", "%s/%s_muscle.fasta" % (alDir, each.split("/")[-1][:-6]), "-quiet"])
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

def kmers(seq,k):
    out=[]
    for i in range(len(seq)-k+1):
        out.append(seq[i:i+k])
    return set(out)

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

def kmerOvlp(one, two, k):
    return len(kmers(one,k).intersection(kmers(two,k)))
    

def pairwise(eD, k, thresh):
    for a, b in it.combinations(eD.keys(), 2):
        if kmerOvlp(eD[a], eD[b], k) >= thresh:
            print("%s\t%s\n%s\t%s\n" % (eD[a], a, eD[b], b))
    
def makeClusts(eD, k, thresh):
    clusts = []
    for n,s in eD.items():
        placed=0
        for i,c in enumerate(clusts):
            matched=0
            for each in list(c.values()):
                if kmerOvlp(s, each, k)>=thresh:
                    clusts[i][n]=s
                    matched=1
            if matched:
                placed+=1
        if not placed:
            clusts.append({n:s})
#        elif placed>1:
#            print (n, placed)
    return(clusts)


def mergeClusts(clusts, k, indivThresh, clustThresh):
    pairs2merge = []
    for i in range(len(clusts)):
        for j in range(i+1, len(clusts)):
            hits=0
            comps=0
            for n1,s1 in clusts[i].items():
                for n2,s2 in clusts[j].items():
                    comps+=1
                    if kmerOvlp(s1, s2, k)>=indivThresh:
                        hits+=1
            if hits/comps >= clustThresh:
                pairs2merge.append([i,j])
    
    groups2merge = combPairs(pairs2merge)
    
    newClusts = []
    singles = []
    merged = []
    for each in groups2merge:
        newClusts.append({})
        for a in each:
            merged.append(a)
            for k,v in clusts[a].items():
                newClusts[-1][k] = v
    for i,info in enumerate(clusts):
        if i not in merged:
            if len(info) == 1:
                singles.append(info)
            else:
                newClusts.append(info)
    return newClusts, singles

def combPairs(pairs):
    d={}
    groups = pairs[::]
    for each in pairs:
        for x in each:
            d[x] = d.get(x, 0) + 1
    for each, count in d.items():
        if count>1:
            groups = combine(each, groups)
    return groups
            
def combine(focal, gList):
    newL = []
    combo = []
    for each in gList:
        if focal in each:
            combo+=each
        else:
            newL.append(each)
    newL.append(list(set(combo)))
    return newL


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

###------------------------------------->>>>    

if __name__ == "__main__":
    main()
