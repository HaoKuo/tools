#!/usr/bin/env python
from __future__ import division

'''
    usage: python filter_genes.py genelist_file target_file output_file 
'''
import os, sys, re, gzip

def main():
    try:
        genelist_file = sys.argv[1]
        target_file = sys.argv[2]
        output_file = sys.argv[3]
        if os.path.exists(target_file) and os.path.exists(genelist_file):
                f1 = open(target_file,"rb")
                header=f1.readline()
                f1.close()
                content = []
                gene_list =[]
                with open(genelist_file,"rb") as f2:
                    for l in f2:
                        gene_list.append(l.strip())
                if os.path.splitext(target_file)[1]=='.gz':
                    for gene in gene_list:
                        mutation_loci_list = grep_gz_file(gene,target_file)
                        content.append(mutation_loci_list)
                else:
                    for gene in gene_list:
                        mutation_loci_list = grep_file(gene,target_file)
                        content.append(mutation_loci_list)
                        print "---Processing {}---".format(gene)
                filtered_file = open(output_file,"wb")
                filtered_file.write(header)
                for c in content:
                    filtered_file.write(''.join(c))
                filtered_file.close()
    except (StandardError),e:
            print("check usage: python filter_genes.py genelist_file target_file output_file")
            sys.exit(0)


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


def grep_file(string, path):
    if os.path.exists(path):
        lines = []
        with open(path,'r') as f:
            for line in f:
                if re.search(string,line):
                    lines +=[line]
    else:
        print('The path \'{}\' is not exist!'.format(path))
    return lines

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


if __name__=='__main__':
    main()


