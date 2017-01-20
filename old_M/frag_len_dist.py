#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
''' Usage: python frag_len_dist.py fragment/ output_name
           files in fragment directory can be plain file or .gz file.
'''
import os, sys, glob
from timeFunc import timeit
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


@timeit
def main():
    args_t = sys.argv[1]
    output_name = sys.argv[2]
    wd = os.getcwd()
    tmp_dir=wd+"/"+args_t
    gzfileList = os.listdir(tmp_dir)
    nameList = [x.split(".")[0] for x in gzfileList]
    dc_list = []
    if os.path.splitext(gzfileList[0])[1] is "gz":
        for gzf in gzfileList:
            dc = pd.read_table(tmp_dir+gzf,compression = 'gzip',header=None,sep='\t',dtype='int32')
            dcol = dc.ix[:,0]
            dcol.columns = gzf.split(".")[0]
            dc_list.append(dcol)
    else:
        for gzf in gzfileList:
            dc = pd.read_table(tmp_dir+gzf, sep='\t', header=None, dtype='int32')
            dcol = dc.ix[:,0]
            dcol.columns = gzf.split(".")[0]
            dc_list.append(dcol)
    combined = pd.concat(dc_list,axis=1,keys=nameList)
    #combined.to_csv("combin_output.txt",index=None,sep="\t")
    df = pd.melt(combined)
    dff = df[df.value<=500]
    g = sns.FacetGrid(dff, col="variable", col_wrap=6)
    g.map(sns.distplot,"value",color = "m")
    g.savefig(output_name+"-insert_length_distribution.jpg")    

if __name__=="__main__":
    main()


