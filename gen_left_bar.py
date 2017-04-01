#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''Usage: python gen_lef_bar.py txt_dir nRows'''

import os, glob, sys, itertools, gzip
from concurrent import futures
import subprocess, multiprocessing
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import math

plt.switch_backend('agg')



def gen_EnrichBar((filename,nRows)):
    print filename, nRows
    data = pd.read_table(filename, sep='\t', nrows=int(nRows))
    #print data.head()
    try :
        data['Term'] = data['Term'].str.replace('^.*~','')
    except :
        data['Term'] = data['Term'].str.replace('^.*:','')
    data['PValue'] = data['PValue'].apply(lambda x: -math.log10(x))

    #print data.head()
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 5))
    sns.set_color_codes("muted")
    ax = sns.barplot(y="Term", x="PValue", data=data,label="Enrichment Score", color="b")
    ax.legend(ncol=2, loc="lower right", frameon=True)
    ax.set(ylabel="", xlabel="-log10(PValue)")
    #plt.yticks(rotation=45)
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    output = os.getcwd() + "/" + os.path.basename(filename).split(".")[0] + "-barplot.png"
    plt.savefig(output, dpi=300)
    
def main():
    txt_dir = sys.argv[1]
    nRows = sys.argv[2]
 #   pool = multiprocessing.Pool(9) 
 #   try :
    para = ((txtFile, nRows) for txtFile in glob.glob(txt_dir+"*.txt"))
    for i, j in para:
        print i,j
        gen_EnrichBar((i,j))
    #pool.map(gen_EnrichBar, para)
 #   except (StandardError), e:
  #      print "Error, no TXT result file found!"
  #      sys.exit(0)

    #with futures.ProcessPoolExecutor(9) as pool:
#        try :
#            
#            para = ((txtFile, nRows) for txtFile in glob.glob(txt_dir+"*.txt"))
#            pool.map(gen_EnrichBar, para)
#        except (StandardError), e:
#            print "Error, no TXT result file found!"
#            sys.exit(0)


if __name__ == "__main__":
    main()

