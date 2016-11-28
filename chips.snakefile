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

rule target:
    input: 
        expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"]), 
        expand("analysis/peaks/{run}/{run}_peaks.bed", run=config['runs'].keys()),
        expand("analysis/fastqc/{sample}_perSeqGC.txt", sample=config["samples"]),
        expand("analysis/peaks/{run}/{run}_sorted_5k_summits.bed", run=config["runs"].keys()),
        expand("analysis/conserv/{run}/{run}_conserv.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_summary.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_stats.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_stats.txt", run=config["runs"].keys()),
        expand("analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam", sample=config["samples"]),
        expand("analysis/peaks/{run}/{run}_4M_peaks.narrowPeak", run=config["runs"].keys()),
        expand("analysis/frips/{run}/{run}_frip.txt",run=config["runs"].keys()),


    message: "Compiling all output"

include: "./modules/align.snakefile"         # rules specific to BWA
include: "./modules/peaks.snakefile"         # peak calling rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
include: "./modules/conservation.snakefile"  # generate conservation plot
include: "./modules/ceas.snakefile"          # annotate peak regions
include: "./modules/frips.snakefile"         # fraction of reads in peaks
include: "./modules/motif.snakefile"         # motif module
include: "./modules/report.snakefile"        # report module
