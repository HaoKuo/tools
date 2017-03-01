'''Supporting functions for the plot of uniformity'''
from __future__ import division

import os, glob, sys, itertools, gzip
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, \
         sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)

def gen_line(detail_txt,filename):
    sns.set_style("darkgrid") 
    plt.switch_backend('agg')
    data = pd.read_table(detail_txt,compression = 'gzip',sep='\t')
    data.columns=['chrom','pos','id','cumPos','depth']
    data2 = DplyFrame(data) >> sift(X["cumPos"] !=0 ) >> select(X.id,X.cumPos,X.depth)
    fig = plt.figure(figsize=(30,12),dpi=600)
    plt.plot(data2['cumPos'],data2['depth'],color="#338844",linewidth=1)
    plt.xlabel('Target position')
    plt.ylabel('Depth')
    plt.title('Coverage across all genomic regions in the panel')
    fig.savefig(filename)


def gen_dot(detail_txt, filename):
    sns.set_style("darkgrid")
    plt.switch_backend('agg')
    data = pd.read_table(detail_txt,compression = 'gzip',sep='\t')
    data.columns=['chrom','pos','id','cumPos','depth']
    data2 = DplyFrame(data) >> sift(X["cumPos"] !=0 ) >> select(X.id,X.cumPos,X.depth)
    fig = plt.figure(figsize=(30,12),dpi=600)
    plt.scatter(data2['cumPos'],data2['depth'],color="#338844",alpha=.4)
    plt.xlabel('Target position')
    plt.ylabel('Depth')
    plt.title('Coverage across all genomic regions in the panel')
    fig.savefig(filename)


