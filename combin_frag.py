#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
''' Usage: python combin_frag.py fragment/ mapped/ output_name
           files in fragment directory can be plain text file or .gz file.
'''
import os, sys, glob
from timeFunc import timeit
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pysam
import itertools
import multiprocessing



def main():
    args_t = sys.argv[1]
    args_m = sys.argv[2]
    output_name = sys.argv[3]

    wd = os.getcwd()
    frag_dir=wd+"/"+args_t
    mapped_dir = wd +"/"+args_m
    rc=gen_fragment_gz(mapped_dir,frag_dir)
    frag(frag_dir, output_name)


def gen_fragment_gz(mapped_dir,frag_dir):
    bam_list = glob.glob(mapped_dir+"*.bam")
    os.chdir(frag_dir)
    #print bam_list
    pool = multiprocessing.Pool(9)
    #pool.map(gen_frag, list(bam_list))
    status = []
    for b in bam_list:
        result = pool.apply_async(gen_frag,(b,))
        status.append(result)
    pool.close()
    pool.join()
    return status

def gen_frag(b):
    f = gzip.open(os.path.basename(b)+".frag.gz","wb")
    with pysam.AlignmentFile(b,"rb") as bam:
        for reads in itertools.islice(bam,0,2000000,2):
            f.write(str(abs(reads.template_length))+'\n')


def frag(tmp_dir,output_name):
    #args_t = sys.argv[1]
    #output_name = sys.argv[2]
    #wd = os.getcwd()
    #tmp_dir=wd+"/"+args_t
    gzfileList = os.listdir(tmp_dir)
    nameList = [x.split(".")[0] for x in gzfileList]
    dc_list = []
    #print os.path.splitext(gzfileList[0])
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


