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
    p.add_option('-m', '--method', default='pearson', help='Correlation method to use. Options: ‘pearson’, ‘kendall’, ‘spearman’ [pearson]')
    p.add_option('-c', '--cmap',  default="viridis", help='Matplotlib colormap to use for plot. [viridis]')
    p.add_option('--dpi', type='int',  default=300, help='Dots per inch to use for plot. [300]')

    opts, args = p.parse_args()
    
    mt.corrMat(opts.input, opts.output, cMap=opts.cmap, setdpi=opts.dpi, meth=opts.method)

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

