#!/usr/bin/env python
"""Len Taing 2023 (TGBTG)
Given a {sample}.virusseq.counts.txt file- cols: chr, start, end, virus, count
this script will:
1. Filter out any rows where the count col = 0
2. for rows that are not 0, the TPM is calculated and outputted

OUTPUT: cols: virus, (raw) count, TPM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [virusseq.count.txt] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--count_file", help="virusseq.count.txt file")
    optparser.add_option("-m", "--mapping_file", help="mapping.txt file that contains samples total read count")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.count_file or not options.mapping_file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #PARSE total number of reads, first line, first col
    total_reads = 1
    with open(options.mapping_file) as f:
        l = f.readline().strip().split(" ")
        total_reads = int(l[0]) / 1000000.0 #SCALING FACTOR = TotalReads per Million

    ls = []
    with open(options.count_file) as f:
        for l in f:
            tmp = l.strip().split("\t")
            count = int(tmp[4])
            #print(count)
            if count > 0:
                virus = tmp[3]
                length_perKb = (int(tmp[2]) - int(tmp[1])) / 1000.0
                rpk = float(count) / length_perKb
                #Add it to the running list
                ls.append((virus, count, rpk))

    #total_rpk = sum([x[2] for x in ls])

    #normalize the RPK values by SCALING FACTOR to get TPM;
    #also sort in desc by TPM
    #ls = sorted([(x[0], x[1], x[2]/total_rpk) for x in ls], key=lambda x: x[2], reverse=True)
    ls = sorted([(x[0], x[1], x[2]/total_reads) for x in ls], key=lambda x: x[2], reverse=True)

    with open(options.output, "w") as out:
        out.write("%s\n" % "\t".join(["Virus","Counts","TPM"]))
        for (v,c,t) in ls:
            out.write("%s\n" % "\t".join([v, str(c), str(t)]))
    
if __name__=='__main__':
    main()


