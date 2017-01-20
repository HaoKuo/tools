#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Parse mapping information from flagstat files.
    usage: python parseMapInf.py flagstat_txt_dir/

'''
from __future__ import division
__author__ = 'Hao Guo'

import os, sys, re, getopt, time
def main():
    flagstatFileDic=sys.argv[1]
    cwd = flagstatFileDic
    os.path.splitext(os.listdir(cwd)[9])
    output = open("all_flagstat.txt", 'w')
    output.write('sample\ttotal_reads\tmapped_reads\tmapped_percent\tproperly_paired\tsingleton_percent\tinterchrom_percent\n')
    for i in os.listdir(cwd):
        a = os.path.splitext(i)
        if len(a) > 1:
            if a[1] == '.txt' and os.path.splitext(a[0])[1] == '.flagstat':
                print "Processing file: %s" %(i)
                sample = re.split(r'\.', i)[0]
                file_path=flagstatFileDic+'/'+i
                flag_file = open(file_path, 'r')
                lines = flag_file.readlines()
                total_reads = re.split(r'\s', lines[0])[0]
                mapped_reads = re.split(r'\s', lines[4])[0]
                mapped_percent = re.search(r'\d+\.\d+', lines[4]).group()
                properly_paired = re.search(r'\d+\.\d+', lines[8]).group()
                singleton_percent = re.search(r'\d+\.\d+', lines[10]).group()
                interchrom_reads = re.split(r'\s', lines[11])[0]
                interchrom_percent = int(interchrom_reads) / int(total_reads) * 100
                #print interchrom_reads, interchrom_percent
                output.write('%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n' % (
                    sample, total_reads, mapped_reads, mapped_percent, properly_paired, singleton_percent,
                    interchrom_percent))
                flag_file.close()
    output.close()


if __name__ == '__main__':
    print '--------------Generating all_flagstat.txt -----------------'
    main()
    print 'Content: \n'
    cat_file = 'cat all_flagstat.txt'
    os.system(cat_file)
    print '-----------------------------------------------------------'




