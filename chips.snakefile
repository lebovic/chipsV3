#!/usr/bin/env python

import os
import sys
import subprocess
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

def addPy2Paths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')

    if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
        config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')


#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = getRuns(config)
addPy2Paths_Config(config)
#-----------------------------------------

def all_targets(wildcards):
    _qdnaseq = config["cnv_qdnaseq_analysis"]
    ls = []
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(peaks_targets(wildcards))
    ls.extend(fastqc_targets(wildcards))
    ls.extend(conservation_targets(wildcards))
    ls.extend(ceas_targets(wildcards))
    ls.extend(frips_targets(wildcards))
    ls.extend(motif_targets(wildcards))
    if _qdnaseq:
        ls.extend(qdnaseq_targets(wildcards))
    return ls

rule target:
    input: 
        all_targets,
        "report.html",

    message: "Compiling all output"
if config['aligner'] == 'bwa':
    include: "./modules/align_bwa.snakefile"     # rules specific to BWA
else:
    include: "./modules/align_bwt2.snakefile"     # rules specific to Bowtie2

include: "./modules/align_common.snakefile"  # common align rules
include: "./modules/peaks.snakefile"         # peak calling rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
include: "./modules/conservation.snakefile"  # generate conservation plot
include: "./modules/ceas.snakefile"          # annotate peak regions
include: "./modules/frips.snakefile"         # fraction of reads in peaks
include: "./modules/motif.snakefile"         # motif module
include: "./modules/contamination.snakefile" # contamination panel module
include: "./modules/qdnaseq.snakefile"       # qdnaseq (CNV) module
include: "./modules/report.snakefile"        # report module
