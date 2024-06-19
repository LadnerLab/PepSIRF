#!/usr/bin/env python

import argparse, os, subprocess
import pandas as pd
# pd.options.mode.copy_on_write = True
# import numpy as np
# import inout as io

# This script reads in a file with NCBI nucleotide accessions and finds any associated Biosample accessions

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-d', '--data',  help='Tab-delimited data matrix containing a column with NCBI nucleotide accession numbers.', required=True)
    reqArgs.add_argument('--outDir',  help='Directory where outputs should be written. One file will be generated per input accession.', required=True)
    reqArgs.add_argument('--outName',  help='Base filename to write Biosample IDs, linked to the nt accessions.', required=True)

    parser.add_argument('--ntName', default="NCBIaccession-NT",  help='Header for input file column with NCBI nucleotide accession numbers.')
    parser.add_argument('--errorLog', default="biosample_errorlog.tsv", help="Name for output file with info on errors encountered.")
    parser.add_argument('--emptyLog', default="biosample_emptyresult.txt", help="Name for output file with info on empty results, likely indicating a lack of associated Biosample.")
    parser.add_argument('--noLog', default="biosample_noaccession.txt", help="Name for output file with info on results without accession annotation.")

    parser.add_argument('rest', nargs=argparse.REMAINDER)
    
    opts = parser.parse_args()

    # Read in input file as a pandas dataframe
    df = pd.read_csv(opts.data, sep="\t", header=0, keep_default_na=False)
    
    if opts.ntName in df:
        accS = df[opts.ntName]
    else:
        print(f"{opts.ntName} was not found in {opts.data}")
        return False

    # Check whether the output directory already exists
    if not os.path.isdir(opts.outDir):
        os.mkdir(opts.outDir)

    #Open a file to use as an error log
    ferr = open(f"{opts.outDir}/{opts.errorLog}", "w")
    ferr.write("Command\tReturnCode\n")
    
    #Open a file to keep track of the samples with empty result
    femp = open(f"{opts.outDir}/{opts.emptyLog}", "w")
    femp.write("Command\n")
    
    #Open a file to keep track of the samples without an accession header in xml
    fnoac = open(f"{opts.outDir}/{opts.noLog}", "w")
    fnoac.write("Command\n")
    
    
    #Dictionary to save biosample accessions
    biosampD = {}
    
    for accRaw in accS:
        if accRaw != "":
            outL=[]
            for acc in accRaw.split(","):
                outName = f"{opts.outDir}/{acc}.xml"
            
                cmd = f"elink -db nuccore -target biosample -id {acc} | efetch -mode xml > {outName}"
                spr = subprocess.run(cmd, shell=True, capture_output=True)
                
                # Flag runs with a non-zero return code
                if spr.returncode != 0:
                    ferr.write(f"{cmd}\t{spr.returncode}\n")
                    outL.append("")
                else:
                    #Check to make sure the file isn't empty
                    if os.stat(outName).st_size == 0:
                        femp.write(f"{cmd}\n")
                        outL.append("")
                        os.remove(outName)
                    else:
                        dfx = pd.read_xml(outName)
                        if "accession" in dfx:
                            outL.append(dfx["accession"][0])
    
                        else:
                            fnoac.write(f"{cmd}\n")
                            outL.append("")
    
            biosampD[accRaw] = ",".join(outL)

    ferr.close()
    femp.close()
    fnoac.close()

    with open(f"{opts.outDir}/{opts.outName}", "w") as fout:
        fout.write(f"{opts.ntName}\tNCBI-Biosample\n")
        for k,v in biosampD.items():
            fout.write(f"{k}\t{v}\n")
    

#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

