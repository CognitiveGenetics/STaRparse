#!#!/usr/bin/env python3

import subprocess, sys, os

    #   MAIN FUNCTION
def main(reads, json, output, savename, globdir, reference,min_sd=3,min_diff=3,min_length=50):
    #   CHECK FILE PATHS
    if output[-1] != "/":
        output += "/"
    # reference file is valid
    if not os.path.isfile(reference):
        raise ValueError(f"Provided reference file does not exist :{reference}")

    #   RUN R FILTER SCRIPTS
    try:
        print("###############     EXTRACT SAMPLE REPEAT DATA     ###############")
        cmd = [os.path.join(globdir,"Source/individual_results/individuals.R"),
                               "-d", globdir,
                                "-r", reads, 
                                "-o", output,
                                "-s", savename,
                                "-R", reference,
                                "-S", str(min_sd),
                                "-D", str(min_diff), 
                                "-L", str(min_length),
                                "-j", json,
        ]
        subprocess.check_call(cmd, shell=False)
        print(savename + " repeat files can be found in directory: " + output)
    except Exception as e:
        sys.exit(print(f"Error:	Failed to extract sample-level repeat data :  {repr(e)}"))

