#!/usr/bin/env python

# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/macs/qc.py
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

def json_macs2_rep(options):
    """
    collect replicates macs2 info to json files
    compared to merged one, collect redundant ratio with --keep-dup 1 option
    """
    input={"all_peak_xls": list(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={"samples":[]}

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    parsed = []
    for i in input["all_peak_xls"]:
        if os.path.exists(i): ## in case only broad peaks would break down sometimes, narrowPeak very seldom no peaks
            parsed.append(_peaks_parse(i))

    if all(list(map(os.path.exists, input['all_peak_xls']))):
        for sample, stat in zip(param["samples"], parsed):
            json_dict["stat"][sample] = stat
        json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action="append", help="input xls files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-b", "--samples", help="paramaters: samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_macs2_rep(options)


if __name__ == '__main__':
    main()
