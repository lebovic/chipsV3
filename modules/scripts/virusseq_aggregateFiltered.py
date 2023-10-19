#!/usr/bin/env python
"""Len Taing 2023 (TGBTG)
Given several {sample}.virusseq.counts.filtered.txt file- 
cols: virus, count, tpm

this script will:
1. aggregate all of the rows across all of the files into a single (pd) table
2. OUTPUT a matrix of raw counts: cols- Sample, virus1, virus2, ..., virusN
3. OUTPUT a matrix of TPM

"""

import os
import sys
import pandas as pd
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [virusseq.count.txt] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="virusseq.count.filtered.txt files")
    optparser.add_option("-s", "--samples", action="append", help="sample names corresponding to reach file")
    optparser.add_option("-o", "--output_count", help="raw count matrix output .csv file")
    optparser.add_option("-t", "--output_tpm", help="tpm matrix output .csv file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.samples or not options.output_count or not options.output_tpm:
        optparser.print_help()
        sys.exit(-1)

    #list of dictionaries- each dict: {VIRUS: count/tpm};
    #list ORDER corresponds to sample order
    counts = []
    tpms = []
    virus_list = ['Sample']
    for (i,ffile) in enumerate(options.files):
        #reads in the filtered.txt file
        with open(ffile) as f:
            ct_dict = {'Sample': options.samples[i]}
            tpm_dict = {'Sample': options.samples[i]}
            #read header
            hdr = f.readline().strip().split("\t")
            for l in f:
                tmp = l.strip().split("\t")
                virus = tmp[0]
                ct = int(tmp[1])
                tpm = float(tmp[2])
                if virus not in virus_list:
                    virus_list.append(virus)
                ct_dict[virus] = ct
                tpm_dict[virus] = tpm
            #Add the dictionaries to the list
            counts.append(ct_dict)
            tpms.append(tpm_dict)
            
    #create a dataframe from the series data
    df_counts = pd.DataFrame(counts, columns=virus_list)
    #NaN values should be set to 0
    df_counts = df_counts.fillna(0)
    #sum up total counts
    total_cts = df_counts.sum(axis=0).drop('Sample', axis=0).sort_values(ascending=False)
    #SORT values by column totals
    df_counts = df_counts[total_cts.index.insert(0,'Sample')]
    #print(df_counts)
    with open(options.output_count, "w") as out:
        df_counts.to_csv(out, sep=",", index=False) #don't print row num

    df_tpm = pd.DataFrame(tpms, columns=virus_list)
    #NaN values should be set to 0
    df_tpm = df_tpm.fillna(0)
    #sum up total counts
    total_tpm = df_tpm.sum(axis=0).drop('Sample', axis=0).sort_values(ascending=False)
    #SORT values by column totals
    df_tpm = df_tpm[total_tpm.index.insert(0, 'Sample')]
    #print(df_tpm)
    with open(options.output_tpm, "w") as out:
        #NOTE: use .round method to round to 3-sig digits
        df_tpm.round(3).to_csv(out, sep=",", index=False) #don't print row num

if __name__=='__main__':
    main()


