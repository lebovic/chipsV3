#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import argparse


def getnew(feature_path,full_table_path):
    #filter:  mrna chromosome protein_coding
    # Read in annotation file
    new_feature = pd.read_csv(feature_path,sep='\t',low_memory=False) 
    new_feature = new_feature.loc[new_feature['seq_type'] == 'chromosome']
    protein_codings = new_feature.loc[new_feature['class']=='protein_coding'].loc[:,'symbol'].tolist()
    new_feature = new_feature.loc[(new_feature['# feature']=='mRNA') & (new_feature['symbol'].isin(protein_codings))]
    new_feature.loc[:,'chromosome'] = ['chr'+str(x) for x in new_feature['chromosome'].tolist()]
    new_feature.loc[new_feature.strand == '+', 'TSS'] = new_feature.start - 1
    new_feature.loc[new_feature.strand == '-', 'TSS'] = new_feature.end - 1
    new_feature.loc[:,'TSS'] = new_feature.TSS.astype(int)
    new_feature.loc[:,'coordinate'] = [x[5]+':'+str(x[7])+'-'+str(x[8]) for x in new_feature.values.tolist()]
    new_feature_bed = new_feature[['chromosome','start','end','coordinate','product_accession','strand','symbol','TSS']] 
    # if refGene_path:
    #     refGene = pd.read_csv(refGene_path,sep = '\t',names = ['bin_num','product_accession','chromosome','strand','start','end','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','symbol','cdsStartStat','cdsEndStat','exonFrames'], low_memory=False)
    #     refGene = refGene[refGene['chromosome'].isin(['chr'+str(i) for i in np.arange(1, 23)]+['chrX', 'chrY'])]
    #     refGene.loc[refGene.strand == '+', 'TSS'] = refGene.start
    #     refGene.loc[refGene.strand == '-', 'TSS'] = refGene.end
    #     refGene.TSS = refGene.TSS.astype(int) 
    #     gene_ann2 = new_feature_bed[new_feature_bed.symbol.isin(refGene.symbol.tolist())]
    #     omit_symbol = [x for x in set(refGene.symbol.tolist()) if x not in set(new_feature_bed.symbol.tolist())]
    #     gene_ann3 = new_feature_bed[new_feature_bed.symbol.isin(omit_symbol)]
    #     # Final full table
    #     full_table = pd.concat([gene_ann2,gene_ann3],axis=0)
    # else:

    # gtf is 1-based system, but bed is 0-based
    new_feature_bed.to_csv(full_table_path, index = None, sep = '\t')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare annotation for RP calculation')
    # parser.add_argument("-r", "--refGene", dest="refGene", default= "" , help="refGene file that downloaded from USCS table")
    parser.add_argument("-f", "--feature", dest="feature", required=True, help="Feature table for updating your annotation")    
    parser.add_argument("-o", "--output", dest="output", required=True, help='prefix for the output file')
    
    args = parser.parse_args()

    # refGene_path=args.refGene
    feature_path=args.feature
    full_table_path=args.output

    getnew(feature_path,full_table_path)




