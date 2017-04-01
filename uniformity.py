'''uniformity.py: plot uniformity figure '''
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

__author__ = 'Hao Guo'

import os, glob, sys, itertools, gzip
import multiprocessing
import subprocess
import argparse
from concurrent import futures

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)


from uniplot import (gen_dot, gen_line)
from timeFunc import timeit


def _cmd_single(args):
    """ Plot coverage per base for single covdep.gz file"""
    inputFile = args.input
    outputFile = args.output
    plotType = args.type

    if plotType in ('line','dot'):
        print "Processing {}".format(inputFile)
        if plotType =='line':
            gen_line(inputFile,outputFile)
        else:
            gen_dot(inputFile, outputFile)
    else:
        raise ValueError("Unknown plot type %r!" % plotType)


def _cmd_batch(args):
    """ Plot a batch of coverage figs for a directory including a bundle of covdep.gz files. """
    inputDir = args.input
    outputDir = args.output
    plotType = args.type
    N_CPU = args.processes
    #print args
    print "Processing directory: {}".format(inputDir)
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    if plotType in ('line','dot'):
        with futures.ProcessPoolExecutor(N_CPU) as pool:
            try :
                para = ((gzfile, outputDir+"/"+os.path.basename(gzfile).split(".")[0]+".jpg", plotType) for gzfile in glob.glob(inputDir+"/*covdep.gz"))
                pool.map(_choosePlotType,para)
            except (StandardError),e:
                print 'Error!'
                sys.exit(0)
    else:
        raise ValueError("Unknown plot type %r!" % plotType)

def _choosePlotType((inputFile,outputFile,plotType)):
    print outputFile
    if plotType =='line':
        gen_line(inputFile,outputFile)
    else:
        gen_dot(inputFile, outputFile)

@timeit
def main():
    sns.set_style("darkgrid")
    plt.rcParams['agg.path.chunksize'] = 100000
    plt.ioff()
    #plt.switch_backend('agg')
    parser = argparse.ArgumentParser(description='Uniformity.py, plot tool for coverage uniformity of sequencing data.', epilog='Contact Hao Guo <guo.hao@genecast.com.cn> for help.')
    subparsers = parser.add_subparsers(help='sub-commands, using -h for more info.')
    
    P_single = subparsers.add_parser('singleplot', help=_cmd_single.__doc__)
    P_single.add_argument('-i', '--input', help="Input covdep.gz file.")
    P_single.add_argument('-o', '--output', help="Output jpg file name.")
    P_single.add_argument('-t', '--type', choices=('line','dot'), default='line',help="Plot type (line,dot).[Default: %(default)s]")
    P_single.set_defaults(func=_cmd_single)


    P_batch = subparsers.add_parser('batchplot', help=_cmd_batch.__doc__)
    P_batch.add_argument('-i', '--input', help="Input directory of covdep.gz file.")
    P_batch.add_argument('-o', '--output', help="Output directory of jpg files.")
    P_batch.add_argument('-t', '--type', choices=('line','dot'), default='line', help="Plot type (line,dot).[Default: %(default)s]")
    P_batch.add_argument('-p','--processes', type=int, nargs='?', default = 1, help="Number of subprocesses to render the coverage fig in parallel.")
    P_batch.set_defaults(func=_cmd_batch)

    args = parser.parse_args()
    args.func(args)
    

if __name__ == '__main__':
    main()




