#!/usr/bin/env python
"""
Len Taing 2024 (TGBTG)

Given a list of files, write to file which ones are empty

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
    optparser.add_option("-f", "--files", action="append", help="list of files")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.output,"w")
    for ffile in options.files:
        size = os.path.getsize(ffile)
        if size == 0:
            out.write(f"{ffile}\n")
    out.close()
    
if __name__=='__main__':
    main()


