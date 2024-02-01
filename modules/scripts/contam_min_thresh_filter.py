#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import pandas as pd
import numpy as np

def main():
    usage = "USAGE: %prog -f [contamination.csv] -o [OUTPUT .csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="contamination.csv files")
    optparser.add_option("-t", "--threshold", help="threshold", default="2.0")
    optparser.add_option("-o", "--output", help="output figure path")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    threshold = float(options.threshold)
    df = pd.read_csv(options.file, index_col=0)
    df = df.astype(float)

    #filter out cols s.t. if the max val in the col < threshold, the col is
    #dropped
    #ref: https://stackoverflow.com/questions/72001582/how-to-remove-columns-that-have-all-values-below-a-certain-threshold
    for col in df[1:]:
        if df[col].max() < threshold:
            df = df.drop([col], axis=1)

    #output
    df.to_csv(options.output)

if __name__=='__main__':
    main()
