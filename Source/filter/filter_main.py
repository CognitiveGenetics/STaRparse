#!/home/dannear/.virtualenvs/py365/bin/python3

import subprocess

    #	MAIN FUNCTION
def main(reads, coverage, output, savename, globdir):
    #	CHECK FILE PATHS
    if output[-1] != "/":
        output += "/"

    #	RUN R FILTER SCRIPTS
    try:
        print("###############     FILTERING OUT POOR COVERAGE READS     ###############")
        subprocess.call(globdir + "/Source/filter/filter.R -d " + globdir + " -r " + reads + " -c " + coverage + " -o " + output + "CGG_Repeats_" + savename + "_filtered.csv", shell=True)
        print("Filtered reads file can be found at: " + output + "CGG_Repeats_" + savename + "_filtered.csv")
    except:
        sys.exit(print("Error:	Failed to filter reads based coverage"))
