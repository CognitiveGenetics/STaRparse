#!#!/usr/bin/env python3

import subprocess, sys

    #   MAIN FUNCTION
def main(reads, output, savename, build, globdir):
    #   CHECK FILE PATHS
    if output[-1] != "/":
        output += "/"

    #   RUN R FILTER SCRIPTS
    try:
        print("###############     SUMMARISING REPEAT DATA     ###############")
        subprocess.call(globdir + "/Source/summaries/summaries.R -d "+ globdir +" -r " + reads +  " -o " + output + " -s " + savename + " -b " + build, shell=True)
        print(savename + " summary files can be found in directory: " + output)
    except:
        sys.exit(print("Error:	Failed to summarize repeat data"))

