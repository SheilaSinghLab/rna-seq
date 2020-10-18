#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import re


paths = sys.argv[1:]
print(paths)
path_samps = paths[0]
path_srr = paths[1]
path_files = paths[2]

gsm_2_samp = {}
with open(path_samps,'r') as f:
    lines = f.readlines()
    for line in range(0,len(lines),2):
        gsm_2_samp[lines[line].rstrip()] = lines[line+1].rstrip()

files = []   
with open(path_files,'r') as f:
    for line in f.readlines():
        files.append(line.rstrip().split('_')[0])
        
srr_2_gsm = pd.read_csv(path_srr)

file_dict = {}
for file in files:
    gsm = srr_2_gsm[srr_2_gsm.Run == file]['Sample Name'].values[0]
    samp = 'WG_Singh_' + gsm_2_samp[gsm].replace('_','') + '_L001_R1_001.fastq'
    if samp in file_dict.values():
        samp =samp.replace('L001','L002')
    file_dict[file+'_1.fastq'] = samp
    
from os.path import join
path='/home/nachmao/projects/def-sheila/nachmao/rna-seq/data/raw_fastq/'
for file in file_dict.keys():
    os.rename(join(path,file), join(path,file_dict[file]))

