#!/usr/bin/env python

import argparse, os, subprocess
import pandas as pd
pd.options.mode.copy_on_write = True
# import numpy as np
# import inout as io

# This script reads in a file with NCBI nucleotide accessions and finds any associated Biosample accessions

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-d', '--data',  help='Tab-delimited data matrix containing a column with NCBI protein accession numbers.', required=True)
    reqArgs.add_argument('--outDir',  help='Directory where outputs should be written. One file will be generated per input accession.', required=True)
    reqArgs.add_argument('--outName',  help='Base filename to write Strain info, linked to the nt accessions.', required=True)

    parser.add_argument('--aaName', default="NCBIaccession-AA",  help='Header for input file column with NCBI protein accession numbers.')
    parser.add_argument('--errorLog', default="strain_errorlog.tsv", help="Name for output file with info on errors encountered.")
    parser.add_argument('--emptyLog', default="strain_emptyresult.txt", help="Name for output file with info on empty results, likely indicating a lack of associated Biosample.")
    parser.add_argument('--noLog', default="strain_noaccession.txt", help="Name for output file with info on results without accession annotation.")

    parser.add_argument('rest', nargs=argparse.REMAINDER)
    
    opts = parser.parse_args()

    # Read in input file as a pandas dataframe
    df = pd.read_csv(opts.data, sep="\t", header=0, keep_default_na=False)
    
    if opts.aaName in df:
        accS = df[opts.aaName]
    else:
        print(f"{opts.aaName} was not found in {opts.data}")
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
    
    with open(f"{opts.outDir}/{opts.outName}", "w") as fout:
        fout.write(f"{opts.aaName}\tStrain\tIsolate\tClone\n")
        for accRaw in accS:
            if accRaw != "":
                strL=[]
                isoL=[]
                cloL=[]
                for acc in accRaw.split(","):
                    if "|" not in acc:
                        outName = f"{opts.outDir}/{acc}.xml"
                        
                        if not os.path.isfile(outName):
                            cmd = f"efetch -db protein -id {acc} -format gb -mode xml > {outName}"
                            spr = subprocess.run(cmd, shell=True, capture_output=True)
                            
                            # Flag runs with a non-zero return code
                            if spr.returncode != 0:
                                ferr.write(f"{cmd}\t{spr.returncode}\n")
                                strL.append("")
                                isoL.append("")
                                cloL.append("")
                                continue
    
                        #Check to make sure the file isn't empty
                        if os.stat(outName).st_size == 0:
                            femp.write(f"{cmd}\n")
                            strL.append("")
                            isoL.append("")
                            cloL.append("")
                            os.remove(outName)
                        else:
                            try:
                                dfx = pd.read_xml(outName, xpath='.//GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier')
                                dfx = dfx.set_index("GBQualifier_name")
                                dfx = dfx['GBQualifier_value']
                                if "strain" in dfx:
                                    if type(dfx["strain"]) == type(pd.Series()):
                                        strL.append(",".join([x for x in dfx["strain"].values if type(x)==type("s")]))
                                    elif type(dfx["strain"])==type("s"):
                                        strL.append(dfx["strain"])
                                    else:
                                        strL.append("")
                                else:
                                    strL.append("")

                                if "isolate" in dfx:
                                    if type(dfx["isolate"]) == type(pd.Series()):
                                        isoL.append(",".join([x for x in dfx["isolate"].values if type(x)==type("s")]))
                                    elif type(dfx["isolate"])==type("s"):
                                        isoL.append(dfx["isolate"])
                                    else:
                                        isoL.append("")
                                else:
                                    isoL.append("")

                                if "clone" in dfx:
                                    if type(dfx["clone"]) == type(pd.Series()):
                                        cloL.append(",".join([x for x in dfx["clone"].values if type(x)==type("s")]))
                                    elif type(dfx["clone"])==type("s"):
                                        cloL.append(dfx["clone"])
                                    else:
                                        cloL.append("")
                                else:
                                    cloL.append("")

        
                            except:
                                strL.append("")
                                isoL.append("")
                                cloL.append("")
                    else:
                        strL.append("")
                        isoL.append("")
                        cloL.append("")

    
                strOut = ",".join(strL)
                isoOut = ",".join(isoL)
                cloOut = ",".join(cloL)
                fout.write(f"{accRaw}\t{strOut}\t{isoOut}\t{cloOut}\n")
        


    ferr.close()
    femp.close()
    fnoac.close()

    

#----------------------End of main()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

