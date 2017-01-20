#!/usr/bin/env python
from __future__ import division

__doc__ = """
SYNOPSIS

    pypy pyVarcheck.py -i <input.mpileup.gz> -o <output.gz> -h <help>
                            or
    pypy pyVarcheck.py -o <output.gz> <input.mpileup.gz>

DESCRIPTION

   Tab delimited output file with:

   chromosome, position, ref_base, depth, mutations, mutation fractions, support reads, chain fractions, support reads for each chain

EXAMPLES

    pypy pyVarcheck.py -h
    pypy pyVarcheck.py -i input.mpileup.gz -o output.varcheck.gz
    pypy pyVarcheck.py -o output.varcheck.gz input.mpileup.gz

"""
__author__ = "Hao Guo"
__version__ = '0.1'

import os, sys, time, getopt, re
import gzip


# read .gz file
def read_gz_file(path):
    if os.path.exists(path):
        with gzip.open(path,'r') as pf:
            for line in pf:
                yield line
    else:
        print('The path \'{}\' is not exist!'.format(path))


# write .gz file
def write_gz_file(path,content):
    try:
        with gzip.open(path,'w') as f:
            f.write(content)
    except:
        print("gz file write error!")


#split file into chunks
def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]


# test '+' existence
def test_in(s,valid):
   test = 0
   for x in valid:
      if x+'+' in s:
          test = 1
          break
   if test == 1:
      return 1
   else:
      return 0


# test '-' existence
def test_del(s,valid):
   test = 0
   for x in valid:
      if x+'-' in s:
          test = 1
          break
   if test == 1:
      return 1
   else:
      return 0


#find first '+'
def find_first_in(s):
    pin = s.index('+')
    numb = ""
    i = 0
    while 1:
        if s[pin + 1 + i].isdigit():
            numb += s[pin + 1 + i]
            i += 1
        else:
            break
    take = int(numb)
    total = s[pin - 1:pin + len(numb) + 1 + take]
    return(total,numb)


#find first '-'
def find_first_del(s):
    pin = s.index('-')
    numb = ""
    i = 0
    while 1:
        if s[pin + 1 + i].isdigit():
            numb += s[pin + 1 + i]
            i += 1
        else:
            break
    take = int(numb)
    total = s[pin - 1:pin + len(numb) + 1 + take]
    return(total,numb)


# merge lower and upper key values
def merge_indels(dic):
    keys = list(set([s.upper() for s in dic.keys()]))
    inserts_ig = {}
    for k in keys:
        try:
            m = dic[k]
        except:
            m=0
        try:
            n = dic[k.lower()]
        except:
            n = 0
        inserts_ig.setdefault(k, m + n)
    return(inserts_ig)


# construct indel list
def construct_indelgeno(inserts_ig,inserts):
    scan=[]
    inlist = map(lambda x: [x, inserts_ig[x]], inserts_ig.keys())
    inlist.sort(lambda x, y: cmp(y[1], x[1]))
    for il in inlist:
        inname = il[0]
        incount = il[1]
        try:
            pcount = inserts[inname]
        except:
            pcount = 0
        try:
            ncount = inserts[inname.lower()]
        except:
            ncount = 0
        scan += [[inname, incount, pcount, ncount]]
    return(scan)


# parse information of the 5th column from mpileup file
def parse_info(info_string, ref, depth):
    inserts = {}
    dels = {}
    quals = {'a': 0, 'A': 0, 'c': 0, 'C': 0, 'g': 0, 'G': 0, 't': 0, 'T': 0, '.': 0, ',': 0, '*': 0}
    valid = ['a', 'A', 'c', 'C', 'g', 'G', 't', 'T', '.', ',', '*']
    x = 0
    temp = info_string
    temp = re.sub(r'\^\S', '', temp)
    #temp = temp.replace('^+', '')
    temp = temp.replace('$', '')
    while test_in(temp, valid) == 1:
        (total,numb) = find_first_in(temp)
        si=2+len(numb)
        inserts.setdefault('+' + total[si:], temp.count(total))
        temp = temp.replace(total, '')
    while test_del(temp, valid) == 1:
        (total,numb) = find_first_del(temp)
        si = 2 + len(numb)
        dels.setdefault('-' + total[si:], temp.count(total))
        temp = temp.replace(total, '')
    while 1:
        if x >= len(temp):
            break
        elif temp[x] in valid:
            quals[temp[x]] += 1
        elif temp[x] == "^":
            x += 1
        elif temp[x] == '+' or temp[x] == '-':
            temp1 = ""
            i = 0
            while 1:
                if temp[x + 1 + i].isdigit():
                    temp1 += temp[x + 1 + i]
                    i += 1
                else:
                    break
            x += int(temp1) + len(temp1)
        x += 1
    del_ig = merge_indels(dels)
    inserts_ig = merge_indels(inserts)
    Aa_HQ = quals['A'] + quals['a']
    Tt_HQ = quals['T'] + quals['t']
    Cc_HQ = quals['C'] + quals['c']
    Gg_HQ = quals['G'] + quals['g']
    match_HQ = quals['.'] + quals[',']
    dls = quals['*']
    if len(inserts_ig) > 0 and len(del_ig) > 0:
        scan1 = construct_indelgeno(inserts_ig, inserts)
        scan2 = construct_indelgeno(del_ig, dels)
        scan = scan1 + scan2
    elif len(del_ig) > 0:
        scan = construct_indelgeno(del_ig, dels)
    elif len(inserts_ig) > 0:
        scan = construct_indelgeno(inserts_ig, inserts)
    else:
        scan = []
    scan += [['A', Aa_HQ, quals['A'], quals['a']], \
             ['T', Tt_HQ, quals['T'], quals['t']], \
             ['C', Cc_HQ, quals['C'], quals['c']], \
             ['G', Gg_HQ, quals['G'], quals['g']], \
             ['*', dls, 0, 0]]
    for w in scan:
        if w[0] == ref:
            w[1] += match_HQ
            w[2] += quals['.']
            w[3] += quals[',']
    scan.sort(lambda x, y: cmp(float(y[1]), float(x[1])))
    scan2 = map(lambda x: [x[0], \
                           round(x[1] / float(depth), 4), \
                           round(x[2] / float(x[1]), 4) if x[1] != 0 else 0, \
                           round(x[3] / float(x[1]), 4) if x[1] != 0 else 0], scan)
    scan11 = filter(lambda x: x[1] > 0, scan)
    scan22 = filter(lambda x: x[1] > 0, scan2)
    name = '|'.join([x[0] for x in scan11])
    freq = '|'.join([str(x[1] * 100) for x in scan22])
    support = '|'.join([str(x[1]) for x in scan11])
    chainF = '|'.join([str(x[2] * 100)+','+str(x[3]*100) for x in scan22])
    chainS = '|'.join([str(x[2]) + ',' + str(x[3]) for x in scan11])
    info = name + '\t' + freq + '\t' + support + '\t' + chainF + '\t' + chainS + '\n'
    return info


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
        # print opts
        # print ('args:%s'%(args))
    except (getopt.GetoptError) as e:
        print('Option error, please use -h for help!')
        sys.exit(0)
    input = None
    for op, value in opts:
        if op == "-h":
            print('Usage: pypy pyVarcheck.py -i <input.mpileup.gz> -o <output.gz> -h <help>')
            print(' Options:\n -h\t: help message\n -i\t: input mpileup file\n -o\t: output varcheck file')
            print(__doc__)
            sys.exit(0)
        elif op == "-o":
            output = value
           # print output
        elif op == "-i":
            input = value
        else:
            print('Wrong option specified.')
            sys.exit(0)
    try:
        output_file = output
    except:
        output_file = re.split(r'\.', args[0])[0] + '.varcheck.gz'
    try:
        fileDic = input or args[0]
        # print fileDic
        if fileDic != "":
            con = read_gz_file(fileDic)
            with gzip.open(output_file, 'w') as f:
                if getattr(con, '__iter__', None):
                    i=0
                    for line in con:
                        chrom = re.split(r'\s', line)[0]
                        pos = re.split(r'\s', line)[1]
                        ref = re.split(r'\s', line)[2].upper()
                        depth = re.split(r'\s', line)[3]
                        bases = re.split(r'\s', line)[4]
                        if int(depth)>0:
                            info = parse_info(bases, ref, depth)
                            # print("%s\t%s\t%s\t%s\t%s" % (chrom, pos, ref, depth, info))
                            f.write("%s\t%s\t%s\t%s\t%s" % (chrom, pos, ref, depth, info))
                            i += 1
                        else:
                            f.write("%s\t%s\t%s\t%s\n" % (chrom, pos, ref, depth))
    except IndexError as e:
        print('Error, try to use ""pypy pyVarcheck.py -i /path/file.mpileup.gz ""')
        sys.exit(1)

if __name__ == '__main__':
    starttime = time.time()
    # print('--------------Parsing information from mpileup file  -----------------')
    main()
   # print('%s was done!'%(sys.argv[1]))
    print('Processing time is : %s' % (time.time() - starttime))
