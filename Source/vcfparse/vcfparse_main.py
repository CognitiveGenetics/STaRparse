#!#!/usr/bin/env python3

    #	IMPORT LIBRARIES
import glob, pandas

    #	OUTPUT FUNCTION
def Output(vcf_data, output, savename):
    df = pandas.DataFrame(vcf_data)
    df.to_csv(output+"Repeats_" + savename  + ".csv", sep='\t', index=None, header=True)
    print("Output CSV file can be found at: "+output+"Repeats_" + savename  + ".csv")
    print("###############     COMPLETE     ###############")

    #	MAIN FUNCTION
def main(vcfinput, output, savename):
    #	CHECK FILE PATHS
    if output[-1] != "/":
        output += "/"
    if vcfinput[-1] != "/":
        vcfinput += "/"
    #	FETCH VCF FILES
    files = glob.glob(vcfinput+'**/*.vcf')
    if files == []:
        files = glob.glob(vcfinput+'*.vcf')
    if files == []:
        print("WARNING: \t No VCF files were found at the specified directory.")
        exit(1)

    #	RUN VCF PARSE FUNCTION
    from From_ExpansionHunter import ExpansionHunter
    ExpansionHunter(files, output, savename)
 

