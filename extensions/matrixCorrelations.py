#!/usr/bin/env python

import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import optparse

#My functions
import matrixtools as mt

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-i', '--input',  help='Input matrix file. [None, REQ]')
    p.add_option('-o', '--output', help='Name for output correlation plot [None, REQ]')
    p.add_option('-c', '--cmap',  default="viridis", help='Matplotlib colormap to use for plot. [viridis]')

    opts, args = p.parse_args()
    
    mt.corrMat(opts.input, opts.output, cMap=opts.cmap)

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

