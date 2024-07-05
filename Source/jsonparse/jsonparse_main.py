#!/usr/bin/env python3

#   IMPORT MODULES
import json, pandas as pd, glob, os

#    DEFINE MAIN FUNCTION
def main(input_path, output_path, savename):
    #           ENSURE CORRECT FILE PATH
    if output_path[-1] != "/": output_path += "/"
    if input_path[-1] != "/": input_path += "/"

    #    COLLECT JSON FILES FROM SPECIFIED INPUT PATH
    files = glob.glob(input_path+'**/*.json')
    if files == []:
        files = glob.glob(input_path+'*.json')
    if files == []:
        print("WARNING:   There are no .json files in this directory")
        exit(1)
    #    INITIALIZE STOREAGE FILE
    f = open(os.path.join(output_path,f"ReadCoverage_{savename}.csv"), 'w+')
    f.write("Sample_ID,Call_ID,MaxSpanningRead,MaxFlankingRead,MaxInrepeatRead".replace(',', '\t'))
    f.close()
    count=0

    #    LOOP THROUGH ALL JSON FILES AT GIVEN DIRECTORY
    print("###############     IMPORTING JSON DATA     ###############")
    for f in files:
        with open(f) as json_file:
            data = json.load(json_file)["LocusResults"]
        base = str(os.path.splitext(os.path.basename(f))[0].split('.json')[0])
        count += 1
        print("Processing JSON file: ", base, "\tFile no.: ", count)

        #    LOOP THROUGH EACH READ WITHIN THE CURRENT JSON FILE
        for x in data:
            repeat = data[x]
            #    ENSURE THAT READ WAS DECTECTED AND DATA IS PRESENT
            if len(repeat) < 5:
                continue
            else:
                #    COLLECT THE LENGTH AND NUMBER OF FLANKING, SPANNING, AND INREPEAT READS
                span = str(repeat['Variants'][x]['CountsOfSpanningReads'])
                flank = str(repeat['Variants'][x]['CountsOfFlankingReads'])
                inread = str(repeat['Variants'][x]['CountsOfInrepeatReads'])

                strings = [base, x, span, flank, inread]
                f = open(output_path+"ReadCoverage_"+savename+".csv", 'a+')
                f.write('\n')
                f.write('\t'.join(strings))
    f.close()
