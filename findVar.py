#!/usr/bin/env python
from __future__ import division
'''
    usage: pypy findVar.py -f varlist.txt -v var_dir -o output_file
           pypy findVar.py -c chr -l location -d input_file_dic/ -o output_file
'''
import os, sys, time, getopt, re, glob
import gzip


def grep(s,pattern):
    return '\n'.join(re.findall(r'^.*%s.*?$'%pattern,s,flags=re.M))

# grep .gz file
def grep_gz_file(string, path):
    if os.path.exists(path):
        lines = []
        with gzip.open(path,'r') as pf:
            for line in pf:
                if re.search(string,line):
                    lines +=[line]
    else:
        print('The path \'{}\' is not exist!'.format(path))
    return lines


# grep list
def grep_list(string ,list):
    if len(list):
        lines = []
        for i in list:
            if re.search(string, i):
                lines +=[i]
    return lines


# grep chr+loc
def grep_chr_loc(chr, loc, gzfile):
    if os.path.exists(gzfile):
        lines = []
        with gzip.open(gzfile,'r') as pf:
            for line in pf:
                if re.search(loc,line):
                    lines +=[line]
        if len(lines) and len(chr):
            lines2 = []
            for ll in lines:
                if re.search(chr, ll):
                    lines2 += [ll]
            return lines2
        else:
            lines2=[]
            return lines2
            print('Location not found or no Chr specified!')
    else:
        print('The path \'{}\' is not exist!'.format(gzfile))


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:d:f:l:o:v:")
    except (getopt.GetoptError),e:
        print "Usage: pypy findVar.py -f varlist.txt"
        sys.exit(0)
    fileDic = '/varcheck/'
    work_Path = os.getcwd()
    var_dic = work_Path + fileDic
    input_dic = var_dic
    output_file = 'vars_in_samples.txt'
    loc = ''
    chr = ''
    varlist_file = ''
    output_file = work_Path + "/" + output_file
    for op, value in opts:
        if op == "-c":
            chr = value
        elif op == "-h":
            print "Please use -f varlist.txt to generate the sample-variants info file."
            print '----------------------------------------'
            sys.exit(0)
        elif op == "-l":
            loc = value
        elif op == "-f":
            varlist_file = work_Path + '/' + value
        elif op == "-d":
            input_dic = work_Path+'/'+ value
        elif op == "-v":
            var_dic = work_Path + '/'+value
        elif op == "-o":
            output_file = value
        else:
            print "Sorry, non location site or file were specified."
            print '----------------------------------------'
            sys.exit(0)
    sample_var_dic = {}
    if loc and chr:
        os.chdir(input_dic)
        file_path_name = glob.glob(r'*.gz')
        file_path_name.sort()
        for file1 in file_path_name:
            samplename = file1.split('.')[0]
            print("Looking up : %s "% samplename)
            tmp_vars = grep_chr_loc(chr, loc, file1)
            print tmp_vars
            sample_var_dic.setdefault(file1, tmp_vars)
        with open(output_file, 'w') as f:
            for k in sample_var_dic:
                for value in sample_var_dic[k]:
                    f.write("%s\t%s"% (k,value))
    elif varlist_file:
        input_dic = var_dic
        os.chdir(input_dic)
        file_path_name = glob.glob(r'*.gz')
        file_path_name.sort()
        for file in file_path_name:
            samplename = file.split('.')[0]
            print("Looking up : %s" % samplename)
            with open(varlist_file, 'r') as vf:
                tmp_vars = []
                for l in vf:
                    a = l.rstrip().split('\t')
                    chr = a[0]
                    loc = a[1]
                    tmp_vars += grep_chr_loc(chr, loc, file)
                sample_var_dic.setdefault(samplename, tmp_vars)
        with open(output_file, 'w') as f:
            for k in sample_var_dic:
                for value in sample_var_dic[k]:
                    f.write("%s\t%s" % (k, value))
    else:
        print("give an chr:loc site using -c -l or file using -f")



if __name__=='__main__':
    t0 = time.time()
   # print '----------------------------------------'
    main()
   # print '----------------------------------------'
    print('Processing time is : %s s' % (time.time() - t0))


