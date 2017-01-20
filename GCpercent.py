#!/usr/bin/env python
from __future__ import division
__doc__='''
    usage: python GCpercent.py bam_file bed_file 
'''
import pysam
from sys import argv

print __doc__

bam_file = argv[1]
bed_file = argv[2]

bam = pysam.AlignmentFile(bam_file, 'rb')
with open(bed_file) as bf :
    for line in bf:
        line_parts = line.strip().split()
        chr = line_parts[0]
        start = int(line_parts[1])
        end = int(line_parts[2])
        read_data = bam.fetch(chr,start,end)
        total_bases = 0
        gc_bases = 0
        for read in read_data:
            seq = read.query_sequence
            total_bases += len(seq)
            gc_bases += len([x for x in seq if x in ["G","g","C","c"]])
        if total_bases is 0:
            gc_per = 'No_Reads'
        else:
            gc_per = '{0:.2f}%'.format(float(gc_bases)/total_bases * 100)
        print '{0}\t{1}'.format(line.strip(),gc_per)

