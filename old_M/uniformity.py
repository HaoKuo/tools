'''uniformity.py: plot uniformity figure '''
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

__author__ = 'Hao Guo'

import os, glob, sys, itertools, gzip
import multiprocessing
import subprocess
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)


from uniplot import (gen_dot, gen_line)



def _cmd_single(args):
    """ Plot coverage per base for single covdep.gz file"""
    #print "hello _cmd_single func!"
    #print args.input,args.output,args.type
    inputFile = args.input
    outputFile = args.output
    plotType = args.type
    if plotType in ('line','dot'):
        if plotType =='line':
            gen_line(inputFile,outputFile)
        else:
            gen_dot(inputFile, outputFile)
    else:
        raise ValueError("Unknown plot type %r!" % plotType)



def parse_args(args=None):
    """Parse the command line."""
    return parser.parse_args(args=args)


def main():
    
    sns.set_style("darkgrid")
    plt.rcParams['agg.path.chunksize'] = 20000
    plt.ioff()
    #plt.switch_backend('agg')
    parser = argparse.ArgumentParser(description='Uniformity.py, plot tool for coverage uniformity of sequencing data.', epilog='Contact Hao Guo <guo.hao@genecast.com.cn> for help.')
    subparsers = parser.add_subparsers(help='sub-commands, using -h for more info.')
    
    P_single = subparsers.add_parser('singleplot', help=_cmd_single.__doc__)
    P_single.add_argument('-i', '--input', help="Input covdep.gz file.")
    P_single.add_argument('-o', '--output', help="Output jpg file name.")
    P_single.add_argument('-t', '--type', choices=('line','dot'), default='line',help="Plot type (line,dot).[Default: %(default)s]")
    P_single.set_defaults(func=_cmd_single)

    args = parser.parse_args()
    args.func(args)
    

if __name__ == '__main__':
    main()




