#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 19:46:13 2021

@author: btu
"""

import pandas as pd
import numpy as np
from Bio.Seq import Seq
import shutil
import sys
try:
    barcodeList=np.loadtxt(f'{sys.argv[1]}/bararr.csv',dtype=str,delimiter='\t')[1:,:]
except:
    barcodeList=np.loadtxt(f'{sys.argv[1]}/bararr.csv',dtype=str,delimiter=',')[1:,:]
# barcodesFwd=pd.read_csv('../barcodesFwd.csv',dtype=str,index_col=0)
# barcodesRev=pd.read_csv('../barcodesRev.csv',dtype=str,index_col=0)

# barcodesFwd['barcode']=barcodesFwd.apply(lambda R: Seq(R.Seq[19:27]), axis=1)
# barcodesFwd['barcodeL']=barcodesFwd.apply(lambda R: Seq(R.Seq[19:]), axis=1)
# barcodesRev['barcode']=barcodesRev.apply(lambda R: Seq(R.Seq[19:27]), axis=1)
# barcodesRev['barcodeL']=barcodesRev.apply(lambda R: Seq(R.Seq[19:]).reverse_complement(), axis=1)
# barcodeCounts=np.zeros([len(barcodesRev),len(barcodesFwd)])
# barcodeSeqInfo={}
issingle=0
if len(sys.argv)>2 and int(sys.argv[2])==1:
    issingle=1
fastaList=[]
for j,seq in enumerate(barcodeList):
    
        fastaList.append('>bw_fwd_%.3d\n'%j)
        fastaList.append(str(seq[2]).upper()+'\n')
        if not issingle:
            fastaList.append('>bw_rev_%.3d\n'%j)
            fastaList.append(str(seq[3]).upper()+'\n')
# fastaList.append('>fwd_mask\n')
# fastaList.append(barcodesFwd.Seq[0][:19]+'NNNNNNNN'+barcodesFwd.Seq[0][27:]+'\n')
# fastaList.append('>rev_mask\n')
# fastaList.append(barcodesRev.Seq[0][:19]+'NNNNNNNN'+barcodesRev.Seq[0][27:]+'\n')
with open('/home/ec2-user/ont-guppy/data/barcoding/barcodes_bw.fasta','w') as f:
    f.writelines(fastaList)
with open('/home/ec2-user/ont-new/ont-guppy/data/barcoding/barcodes_bw.fasta','w') as f:
    f.writelines(fastaList)
np.savetxt(f'{sys.argv[1]}/samples.csv',barcodeList[:,0],fmt="%s")
#shutil.copyfile('barcodes_bw.fasta','/home/ec2-user/ont-guppy/data/barcoding/barcodes_bw.fasta')
lastIdx=len(barcodeList[:,0])-1
cfgTemplate=f'''[loading_options]
barcodes_filename =barcodes_bw.fasta 
double_variants_frontrear =true
# ############### PCR barcoding kit ###############

[BTU_%03i]
compatible_kits = bw-4
first_index =0
last_index = {lastIdx}
kit = bw 
normalised_id = bw_%03i
scoring_function = MAX
#mask1 = fwd_mask
#mask2 = rev_mask 
barcode1 = bw_fwd_%03i
barcode2 = bw_rev_%03i'''

cfgTemplatesingle=f'''[loading_options]
barcodes_filename =barcodes_bw.fasta 
#double_variants_frontrear =true
# ############### PCR barcoding kit ###############

[BTU_%03i]
compatible_kits = bw-single
first_index =0
last_index = {lastIdx}
kit = bw 
normalised_id = bw_%03i
scoring_function = MAX
#mask1 = fwd_mask
#mask2 = rev_mask 
barcode1 = bw_fwd_%03i
#barcode2 = bw_rev_%03i'''
if issingle:
    with open ('/home/ec2-user/ont-guppy/data/barcoding/bw-single.cfg','w') as f:
        f.write(cfgTemplatesingle)
else:
    with open ('/home/ec2-user/ont-guppy/data/barcoding/bw-3.cfg','w') as f:
        f.write(cfgTemplate)
    with open ('/home/ec2-user/ont-new/ont-guppy/data/barcoding/bw-3.cfg','w') as f:
        f.write(cfgTemplate)
