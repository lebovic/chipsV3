#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -o [OUTPUT FIGURE]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="sample_mapping.txt files")
    optparser.add_option("-o", "--output", help="output figure path")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    with open(options.file) as ffile:
        #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
        sampleID = options.file.strip().split("/")[-1].split('.')[0]
        #ALSO remove the suffix '_mapping' from the sampleID name
        if sampleID.endswith('_mapping'):
            sampleID = sampleID.replace("_mapping","")
        
        tmp = ffile.readlines()
        total = int(tmp[0].strip().split()[0])
        mapped = int(tmp[6].strip().split()[0])
        uniq_mapped = int(tmp[-1].strip())

        data = pd.Series([total,mapped,uniq_mapped],
                         index=['Total Reads', 'Mapped Reads', 'Uniquely Mapped Reads'])

        sns.set(style="whitegrid",font_scale=1.2)
        f, ax= plt.subplots(figsize = (9, 5))
        ax.set_title('%s Mapped Rate' % sampleID)
        sns.barplot(x=data.values,y=data.index,palette="Blues_d")
        f.savefig(options.output, dpi=100,bbox_inches='tight')       

    

if __name__=='__main__':
    main()

    
