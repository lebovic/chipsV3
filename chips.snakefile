#!/usr/bin/env python

import pandas as pd

def getRuns(config):
    """parse metasheet for Run groupings"""
    #metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')    
    #config['runs'] = metadata.index
    ret = {}
    f = open(config['metasheet'])
    hdr=f.readline().strip().split(',')
    for l in f:
        tmp = l.strip().split(",")
        ret[tmp[0]] = tmp[1:]

    config['runs'] = ret
    return config


#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = getRuns(config)
#-----------------------------------------

rule target:
    input: expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"]), expand("analysis/peaks/{run}/{run}_peaks.bed", run=config['runs'].keys())
    message: "Compiling all output"

include: "./modules/align.snakefile"         # rules specific to BWA
include: "./modules/peaks.snakefile"         # peak calling rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
include: "./modules/conservation.snakefile"  # generate conservation plot
include: "./modules/ceas.snakefile"  # annotate peak regions

