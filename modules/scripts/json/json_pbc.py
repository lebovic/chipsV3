#!/usr/bin/env python

# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/library/qc.py
from optparse import OptionParser
import sys
import os
import json

def json_dump(json_dict):   # json
    """
    dump out uniform json files for collecting statistics
    :param json_dict: output python dict to json
    :return: json_file name
    """
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file

def json_pbc(options):
    """
    collect conservation_plot output Phastcon score
    """
    input={"pbc": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    json_dict = {"stat": {}, "input": input, "output": output}
    
    sampleID = str(os.path.abspath(options.input)).strip().split("/")[-1].split('.')[0]
    #remove the _pbc ending
    sampleID = sampleID.replace("_pbc", "")
    f = open(os.path.abspath(options.input))
    firstLine = f.readline().strip().split()
    #print(firstLine)
    N1 = int(firstLine[1]) #record number of reads with just one location
    #NOW SUM over the rest of the reads
    
    Nd = N1
    for l in f:
        tmp = l.strip().split()
        Nd += int(tmp[1])

    json_dict["stat"][sampleID] = {}
    json_dict["stat"][sampleID]["N1"] = N1
    json_dict["stat"][sampleID]["Nd"] = Nd
    json_dict["stat"][sampleID]["PBC"] = round(N1/Nd, 3)
    json_dict["param"]={"samples":sampleID}
    json_dump(json_dict)

def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    (options, args) = optparser.parse_args(sys.argv)
    json_pbc(options)


if __name__ == '__main__':
    main()
