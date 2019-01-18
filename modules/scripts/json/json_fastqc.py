#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/fastqc/qc.py
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

def json_fastqc(options):
    input={"fastqc_summaries": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={"ids": options.ids,
           "id": options.ID}
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    stat = json_dict["stat"]
    for a_summary, a_id in zip(input["fastqc_summaries"], param["ids"]):
        parsed = _fastqc_parse(input=a_summary)
        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["sequence_length"] = parsed["sequence_length"]
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    optparser.add_option("-s", "--ids", help="sample bases")
    optparser.add_option("-I", "--ID", help="sample ID")
    (options, args) = optparser.parse_args(sys.argv)
    json_fastqc(options)


if __name__ == '__main__':
    main()
