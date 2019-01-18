# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/frip/qc.py

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

def json_frip(options):
    input={"frip": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={"samples":options.samples}
    """
    input is *.frip
    output is conf.json_prefix + "_frip.json"
    param for matching samples
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    for i,s in zip(input["frip"], param["samples"]):
        inf = open(i).read().strip().split(",")
        json_dict["stat"][s] = {}
        json_dict["stat"][s]["info_tag"] = int(inf[0])
        json_dict["stat"][s]["total_tag"] = int(inf[1])
        json_dict["stat"][s]["frip"] = float(int(inf[0]))/int(inf[1])
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action="append", help="input frip files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-s", "--samples", help="paramaters: samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_frip(options)


if __name__ == '__main__':
    main()
