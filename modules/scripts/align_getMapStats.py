#!/usr/bin/env python
"""Script to collect the mapping statistics from across all samples. 
outputs to stdout:
Sample,Mapped,Total,Percentage
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -f [FPKM FILE_2] ...-f [FPKM FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of sample_mapping.txt files (note: these are snakemake temp files)")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    print(",".join(["Sample","Total","Mapped","UniquelyMapped"]))

    for f in options.files:
        #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
        sampleID = f.strip().split("/")[-1].split('.')[0]
        #ALSO remove the suffix '_mapping' from the sampleID name
        if sampleID.endswith('_mapping'):
            sampleID = sampleID.replace("_mapping","")

        with open(f) as ffile:
            #Total reads needs to be parsed from 1st line; mapped is 7th line
            #uniq is last line
            tmp = ffile.readlines()
            total = int(tmp[0].strip().split()[0])
            mapped = int(tmp[6].strip().split()[0])
            uniq_mapped = int(tmp[-1].strip())
            print(",".join([sampleID,str(total),str(mapped),str(uniq_mapped)]))

if __name__=='__main__':
    main()


