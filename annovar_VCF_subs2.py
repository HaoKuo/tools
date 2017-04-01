#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''parse mapping information from flagstat files
'''

from __future__ import division
__author__ = 'Hao Guo'
'''
    __doc__
    Generate annotation (1000 genomes 2015 version) of VCF files.
    e.g.: python annovar_VCF_subs.py vcf_dir
'''
import os, sys, re, getopt, time
import multiprocessing


#os.chdir("E:\\share")
#os.chdir("E:\\share\\annovar")
def main():
    try:
        fileDic =  sys.argv[1]
        if fileDic!="" and os.path.exists(fileDic):
           # print fileDic
            #work_Path = os.getcwd()
            annovar_Path = fileDic
           # print annovar_Path
            os.chdir(annovar_Path)
            cwd = os.getcwd()
            files = os.listdir(cwd)
            files.sort()
            output_dic_path = cwd + '/anno/'
            if not os.path.exists(output_dic_path):
                os.mkdir('anno')
            samples = [re.split(r'\.', f)[0] for f in files]
            file_path_name = [os.path.join(cwd, d) for d in files]
            command1='/work-a/user/zhanggh/Software/ANNOVAR/annovar/table_annovar.pl '
            dbDir=' /work-a/user/zhanggh/Software/ANNOVAR/annovar/humandb'
            tmp_str1=' -buildver hg19 -out '
            tmp_str2=' -remove -protocol refGene,cytoBand,genomicSuperDups,clinvar_20160302,cosmic70,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,ljb26_all,exac03 -operation g,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput'
            pool = multiprocessing.Pool(min(20,len(file_path_name)))

            for f in file_path_name:
                if os.path.isfile(f) and os.path.splitext(f)[1]=='.vcf':
                    str = command1 + f + dbDir + tmp_str1 + output_dic_path
                    s1 = re.split(r'\/', f)[-1]
                    s2 = re.split(r'\.',s1)[-2]
                    s3 = re.split(r'\.', s1)[0]
                    command_line=str + s2 + '_' + s3 + tmp_str2
                    print '-----------------------------------------------------------'
                    print "Processing file: " + f
                    print "Command line is :"
                    print command_line
                    pool.apply_async(os.system,(command_line,))
                    #os.system(command_line)
            pool.close()
            pool.join()
            print "Sub-processes done!"
    except (StandardError),e:
        #logging.exception(e)
        print "Option error, Do use:  ""python annovar_VCF_subs.py /dir/"""
        print '----------------------------------------'
        sys.exit(0)

if __name__ == '__main__':
    print '--------------Running Annovar programme -----------------'
    main()
    #print 'Content: \n'
    # cat_file = 'cat all_flagstat.txt'
    # os.system(cat_file)
    print '-----------------------------------------------------------'
