#!/usr/bin/env python
"""Script to collect the peak statistics from across all samples. 
#Total Peaks, #10FC peaks, #20FC peaks
Sample,Total,10FC,20FC
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of _sorted_peaks.narrowPeak files")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.output,"w")
    out.write("%s\n" % ",".join(["Sample","Total","10FC","20FC"]))

    for f in options.files:
        #TRY to infer the SAMPLE NAMES
        sampleID = f.strip().split("/")[-1].split('.')[0]
        if sampleID.endswith('_sorted_peaks'):
            sampleID = sampleID.replace("_sorted_peaks","")

        f = open(f)
        #start the counts
        tot = fc_10 = fc_20 = 0
        for l in f:
            tmp = l.strip().split("\t") 
            #note FC is 7th col
            fc = float(tmp[6])
            if fc >= 20.0:
                fc_20 += 1
            if fc >= 10.0:
                fc_10 += 1
            tot += 1

        out.write("%s\n" % ",".join([sampleID,str(tot),str(fc_10),str(fc_20)]))
        f.close()

if __name__=='__main__':
    main()


