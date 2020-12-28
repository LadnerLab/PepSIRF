#!/usr/bin/env python

import optparse
import inout as io

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--inp',  help='Input linkage map. [None, REQ]')
    p.add_option('-t', '--taxIDs',  help='File containing a list of taxonomic IDs to remove from the linkage map. [None, REQ]')
    p.add_option('-o', "--out", help='Name for output file. [None, REQ]')

    opts, args = p.parse_args()

    toExclude = io.fileEmptyDict(opts.taxIDs, header=False)

    with open(opts.out, "w") as fout:
        with open(opts.inp, "r") as fin:
            lc=0
            for line in fin:
                lc+=1
                if lc==1:
                    fout.write(line)
                else:
                    cols = line.rstrip("\n").split("\t")
                    if cols[1]:
                        links = [x.split(":") for x in cols[1].split(",")]
                        links = [":".join(x) for x in links if x[0] not in  toExclude]
                        fout.write("%s\t%s\n" % (cols[0], ",".join(links)))
                    else:
                        fout.write(line)

#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()
