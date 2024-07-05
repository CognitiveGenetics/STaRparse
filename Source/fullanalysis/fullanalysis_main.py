#!/usr/bin/env python3

    #	IMPORT LIBRARIES / MODULES
import sys, glob, pandas, json, os
#globdir = os.path.dirname(os.path.realpath(__file__))

    #	MAIN FUNCTION
def main(vcfinput, jsoninput, output, savename, build, humandb, globdir, reference=None):
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
    reads_file_path = os.path.join(output,f"Repeats_{savename}.csv")
    coverage_file_path = os.path.join(output,f"ReadCoverage_{savename}.csv")
    main(reads_file_path, coverage_file_path, output, savename, globdir)

    #	RUN SUMMARY FUNCTIONS
    sys.path.insert(1, globdir + "/Source/summaries")
    from summaries_main import main
    reads_filtered_path = output + "Repeats_" + savename + "_filtered.csv"
    try:
        main(reads_filtered_path, output, savename, build, humandb, globdir)
    except Exception as e:
        print(f"WARNING: Could not build summaries: {repr(e)}")
        if not reference:
            print("ERROR : Cannot continue without reference file")
            sys.exit(1)
        pass

    #	EXTRACT SAMPLE FUNCTIONS
    sys.path.insert(1, globdir + "/Source/individual_results")
    from individuals_main import main
    if not reference:
        reference = os.path.join(output, f"{savename}_by_locus.csv")
    reads_filtered_path = output + "Repeats_" + savename + "_filtered.csv"
    main(reads_filtered_path, jsoninput, output, savename, globdir, reference)


