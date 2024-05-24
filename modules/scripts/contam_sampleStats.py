#!/usr/bin/env python
"""
Len Taing 2024 (TGBTG)

Script to collect the contamination statistics from across all species. 
The scripts two paralel lists: a list of contamination stat files and their
associated species name

Output:
{species0} {contents of file i.e. % contam}
{species1} {contents of file i.e. % contam}
...
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -n [species 1] ... -n [species 2] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of contamination stats files")
    optparser.add_option("-n", "--names", action="append", help="list of associated species name")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.names or not options.output:
        optparser.print_help()
        sys.exit(-1)

    if len(options.files) != len(options.names):
        print("ERROR: the number of files must match the number of names!")
        sys.exit(-1)

    out = open(options.output,"w")
    for (ffile, name) in zip(options.files, options.names):
        with open(ffile) as f:
            percent = f.read().strip()
            out.write(f"{name} {percent}\n")

if __name__=='__main__':
    main()


