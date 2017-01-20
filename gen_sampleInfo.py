#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' Generate sample info file for covstat.cpp'''

__author__='Hao Guo'

import os,sys,re,getopt,time

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:o:")
    except (getopt.GetoptError),e:
        print "Option error, please use -h for help!"
        print '----------------------------------------'
        sys.exit(0)

    panel_name=''
    output_file='sample_info.txt'
    for op, value in opts:
        if op == "-p":
            panel_name = value
        elif op == "-h":
            print "Please use -p panel_name to generate the sample info file."
            print '----------------------------------------'
            sys.exit(0)
        elif op == "-o":
            output_file=value
        else:
            print "Sorry, non panel names were specified."
            print '----------------------------------------'
            sys.exit(0)

    work_Path=os.getcwd()
    mpileup_Path=work_Path+'/mpileup_all'
    os.chdir(mpileup_Path)

    cwd = os.getcwd()
    files = os.listdir(cwd)
    files.sort()
    samples = [re.split(r'\.', f)[0] for f in files]
    filePathName=[os.path.join(cwd, d) for d in files]
 
    if len(panel_name)>0:
        print "Input panel:",panel_name
    else:
        print "No panel was given."
    os.chdir(work_Path)
    sample_info_file=open(output_file,'w')
    for i in range(len(samples)):
        sample_info_file.write('%s\t%s\t%s\n' % (samples[i],panel_name,filePathName[i]))
    sample_info_file.close()
    print "Generating sample info file: %s "%(output_file)
    cat_file='cat '+output_file
    os.system(cat_file)

if __name__=='__main__':
    t0 = time.clock()
    print '----------------------------------------'
    #mpileup_Path=os.getcwd()+'/mpileup_target'
    #os.chdir(mpileup_Path)
    main()
    print '----------------------------------------'
    print 'Processing time is ' + str(time.clock() - t0) + ' seconds'




