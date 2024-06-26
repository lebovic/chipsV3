#!/usr/bin/env python

import os
import sys
import subprocess
from optparse import OptionParser

def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-o", "--output", help="output files")
    (options, args) = optparser.parse_args(sys.argv)
    out = subprocess.Popen("chromap -v", shell=True, stderr=subprocess.PIPE).stderr.read().decode('utf-8')
    version = out.strip()
    with open(options.output,"w+") as output:
        output.write(version)

if __name__=='__main__':
    main()

