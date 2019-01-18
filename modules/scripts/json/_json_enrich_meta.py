#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/enrichment/dc.py
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

def json_enrich_meta(options):
    """ enrichment in meta regions
    """
    input = {'meta':str(os.path.abspath(options.input)), 'mapped':str(os.path.abspath(options.mapped))}
    output = {"json": str(os.path.abspath(options.output))}
    param = {'dhs': str(os.path.abspath(options.dhs)), 'down': option.down, 
             'has_dhs':str(os.path.abspath(options.HasDHS)), 'id':options.ID, 'samples':options.samples}
    json_dict = {"stat": {}, "input": input, "output": output, "param":param}
    for n, s in enumerate(param['samples']):
        ## total mapped reads
        mapped = float(open(input["mapped"][n]).readlines()[2].split()[0])
        json_dict['stat'][s] = {}
        meta = open(input['meta'][n]).read().strip().split(",")
        meta = list(map(float, meta))
        if not param["down"]:
            json_dict['stat'][s]['exon'] = meta[0]/mapped
            json_dict['stat'][s]['promoter'] = meta[1]/mapped ## use all mapped reads
        else:
            json_dict['stat'][s]['exon'] = meta[0]/meta[2]
            json_dict['stat'][s]['promoter'] = meta[1]/meta[2] ## use 4M reads

        if param['has_dhs']:
            dhs = open(param["dhs"][n]).read().strip().split(",")
            dhs = list(map(float, dhs))
            if not param["down"]:
                json_dict['stat'][s]['dhs'] = dhs[0]/mapped
            else:
                json_dict['stat'][s]['dhs'] = dhs[0]/dhs[1]
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input meta files")
    optparser.add_option("-m", "--mapped", help="input mapped files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-d", "--dhs", help="set dhs")
    optparser.add_option("-H", "--HasDHS", help="set has dhs")
    optparser.add_option("-D", "--down", help="set down")
    optparser.add_option("-I", "--ID", help="set ID")
    optparser.add_option("-s", "--samples", help="set samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_enrich_meta(options)


if __name__ == '__main__':
    main()
