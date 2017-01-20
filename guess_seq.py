#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
    usage: python guess_seq2.py
'''
import os, sys, time, re, glob, itertools
import subprocess, multiprocessing, shlex
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import pysam



fastq_dir = "/work-z/user/guoh/tech-RD/molecular_barcode/fastq/"
test_fq1 = "/work-z/user/guoh/tech-RD/molecular_barcode/fastq/S370_05B_CHG015512-Sample20161103-1-I13-AACGTGAT_L007_R1.fastq.gz"


bam_dic = "/work-z/user/guoh/tech-RD/molecular_barcode/mapped/"
bam_file1 = "/work-z/user/guoh/tech-RD/molecular_barcode/mapped/S1_I13.bam"
bam_file2 = "/work-z/user/guoh/tech-RD/molecular_barcode/mapped/S1_I14.bam"

def main():
    guess_barcode(bam_file1,"S1_I13_1M",1000000)
    guess_barcode(bam_file2,"S1_I14_1M",1000000)


def guess_barcode(bam_file,output,Nreads):

    ''' usages: guess_barcode(bam_file1,"S1_I13_1M",1000000)    
'''    
    Frag_dic ={}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        #for reads in bam:
        for reads in itertools.islice(bam,0,Nreads):
            if 1<abs(reads.template_length)<150:
                if reads.template_length>0:
                    Frag_dic.setdefault(reads.query_name,{}).setdefault('L',[reads.cigartuples,reads.query_sequence])
                else:
                    Frag_dic.setdefault(reads.query_name,{}).setdefault('R',[reads.cigartuples,reads.query_sequence])
    #print Frag_dic
    #print len(Frag_dic)
    seq_L_list = []
    seq_R_list = []
    for key in itertools.islice(Frag_dic,0,None):
        if len(Frag_dic[key])>1:
            if Frag_dic[key]['L'][0][0][0]==4:
                tmplen = Frag_dic[key]['L'][0][0][1]
                tmpseq = Frag_dic[key]['L'][1]
                seq_L_list.append(tmpseq[0:tmplen])
            if Frag_dic[key]['R'][0][-1][0]==4:
                tmplen = Frag_dic[key]['R'][0][-1][1]
                tmpseq = Frag_dic[key]['R'][1]
                seq_R_list.append(tmpseq[-tmplen:])
    #print seq_L_list
    #print seq_R_list
    La= max([len(i) for i in seq_L_list])
    Ra= max([len(i) for i in seq_R_list])
    minLR = min(La,Ra)
    NminLR = -minLR-1
    #print NminLR
    LL=[x.rjust(max([len(i) for i in seq_L_list])) for x in seq_L_list]
    RR=[x.ljust(max([len(i) for i in seq_R_list])) for x in seq_R_list]
    R_seq_obl={}
    L_seq_obl={}
    for po in range(minLR):
        str1 = reduce(lambda x,y:x+y,map(lambda x:x[po],RR))
        dict_char_tmp = {s:str1.count(s) for s in str1}
        ordered_baseFreq_list = sorted(dict_char_tmp.items(),key=lambda item:item[1],reverse=True)
        R_seq_obl.setdefault(po,ordered_baseFreq_list)
    for po in range(-1,NminLR,-1):
        str2 = reduce(lambda x,y:x+y,map(lambda x:x[po],LL))
        dict_char_tmp = {s:str2.count(s) for s in str2}
        ordered_baseFreq_list = sorted(dict_char_tmp.items(),key=lambda item:item[1],reverse=True)
        L_seq_obl.setdefault(po,ordered_baseFreq_list)
    #print L_seq_obl
    #print R_seq_obl
    f=open(output+"-ligated_barcode.txt","wb")
    L_seq=[]
    for i in range(-minLR,0,1):
        f.write('{} left base is:  '.format(abs(i)))
        f.write(''.join([x.__str__() for x in L_seq_obl[i]])+"\n")
        if L_seq_obl[i][0][0]!=" ":
            L_seq.append(L_seq_obl[i][0][0])
        else:
            L_seq.append(L_seq_obl[i][1][0])
    #f.write("Left barcode is : " + ''.join(L_seq)+"\n")
    f.write("soft clip sequence number is: {} \n".format(len(seq_L_list)))
    f.write("Left max soft clip length is: {} \n".format(max([len(x) for x in seq_L_list])))
    f.write("Left barcode is : " + ''.join(L_seq)+"\n")
    R_seq=[]
    for i in range(minLR):
        f.write('{} right base is: '.format(i+1))
        f.write(''.join([x.__str__() for x in R_seq_obl[i]])+"\n")
        if R_seq_obl[i][0][0]!=" ":
            R_seq.append(R_seq_obl[i][0][0])
        else:
            R_seq.append(R_seq_obl[i][1][0])
    #f.write("Right barcode is : " + ''.join(R_seq)+"\n")
    f.write("soft clip sequence number is : {} \n".format(len(seq_L_list)))
    f.write("Right max soft clip length is: {} \n".format(max([len(x) for x in seq_R_list])))
    f.write("Right barcode is : " + ''.join(R_seq)+"\n")


if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Processing time is : %s s' % (time.time() - t0))

