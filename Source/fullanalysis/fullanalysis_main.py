#!/home/dannear/.virtualenvs/py365/bin/python3

    #	IMPORT LIBRARIES / MODULES
import sys, glob, pandas, json, os
#globdir = os.path.dirname(os.path.realpath(__file__))

    #	MAIN FUNCTION
def main(vcfinput, jsoninput, output, savename, build, globdir):
    #	CHECK FILE PATHS ARE COMPLETE
    if output[-1] != "/":
        output += "/"
    if vcfinput[-1] != "/":
        vcfinput += "/"
    if jsoninput[-1] != "/":
        jsoninput += "/"

    #	PARSE VCF FILES
    sys.path.insert(1, globdir + "/Source/vcfparse")
    from vcfparse_main import main
    main(vcfinput, output, savename)

    #	PARSE JSON FILES
    sys.path.insert(1, globdir + "/Source/jsonparse")
    from jsonparse_main import main
    main(jsoninput, output, savename)

    # 	FILTER READS
    sys.path.insert(1, globdir + "/Source/filter")
    from filter_main import main
    reads_file_path = output + "CGG_Repeats_" + savename + ".csv"
    coverage_file_path = output + "ReadCoverage_" + savename + ".csv"
    main(reads_file_path, coverage_file_path, output, savename)

    #	RUN SUMMARY FUNCTIONS
    sys.path.insert(1, globdir + "/Source/summaries")
    from summaries_main import main
    reads_filtered_path = output + "CGG_Repeats_" + savename + "_filtered.csv"
    main(reads_filtered_path, output, savename, build)

