#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:16:10 2019

@author: dannear
"""
import os, vcf

def ExpansionHunter(vcffiles, output, savename):
   #           IMPORT VCF FILES AND DESIRED EXTRACT DATA
    print("###############     IMPORTING VCF DATA     ###############")

    count = 0
    # output is directory. 
    try:
       os.makedirs(output,exist_ok=True)
    except Exception as e:
       print(f"Could not make output directory : {repr(e)}")
       raise e
    f = open(os.path.join(output,f"Repeats_{savename}.csv"),"w+")
    f.write('Call,Sample_ID,Chr,Start,End,GT,Ref_Units,Allele1_Units,Allele2_Units'.replace(',', '\t'))
    f.close()

    for x in vcffiles:
        f = open(os.path.join(output,f"Repeats_{savename}.csv"),"a+")
        count += 1
        base = os.path.splitext(os.path.basename(x))[0]
        print("Processing VCF file: ", base, "\tFile no.: ", count)
        my_vcf = vcf.Reader(filename=x)
        for record in my_vcf:
            samp = record.samples
            for i in samp:
                 if i['GT'] == ".":
                     continue
                 Call = record.INFO['REPID']
                 chrom = record.CHROM
                 start = int(record.POS)
                 end = int(record.INFO["END"])
                 units = i['REPCN']
                 gt = i['GT']
                 ref = round(record.INFO["REF"])

                 if '/' not in units:
                     Allele1_Units=units
                     Allele2_Units=0
                 else:
                     Allele1_Units=units.split("/")[0]
                     Allele2_Units=units.split("/")[1]

                 strings = [str(Call), str(base), str(chrom), str(start), str(end), str(gt), str(ref), str(Allele1_Units), str(Allele2_Units)]

                 f.write('\n')
                 f.write("\t".join(strings))

        f.close()

