#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/bwa/qc.py
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

def json_map(options):
    """
    input samtools flagstat standard output
    output json files
    kwargs for matching replicates order
    keep one value for each json for easier loading to html/pdf template
    example:
    3815725 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 duplicates
    3815723 + 0 mapped (100.00%:-nan%)
    """
    input={"bwa_mapped": list(os.path.abspath(options.input)), "bwa_total": list(os.path.abspath(options.total))}
    output={"json": str(os.path.abspath(options.output))}
    param={"sample": options.samples}

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for mapped, total, sam in zip(input["bwa_mapped"], input["bwa_total"], param["sample"]):
        inft = open(total, 'rU')
        infm = open(mapped, 'rU')
        json_dict["stat"][sam] = {}
        json_dict["stat"][sam]["mapped"] = int(infm.readlines()[2].split()[0])
        json_dict["stat"][sam]["total"] = int(inft.readlines()[0].strip())
        inft.close()
        infm.close()
    json_dump(json_dict)



def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input bwa mapped files")
    optparser.add_option("-t", "--total", help="input bwa total files")
    optparser.add_option("-o", "--output", help="output files")
    optparser.add_option("-s", "--samples", help="paramaters: samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_map(options)


if __name__ == '__main__':
    main()
