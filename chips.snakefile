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
        #ALIGN_ALL - note KEEP these insync!
        #expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}.sorted.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.sorted.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}.unmapped.fq.gz", sample=config["samples"]),
        "analysis/align/mapping.csv",
        #PEAKS_ALL
        expand("analysis/peaks/{run}/{run}_peaks.bed", run=config['runs'].keys()),
        expand("analysis/peaks/{run}/{run}_sorted_peaks.bed", run=config["runs"].keys()),
        expand("analysis/peaks/{run}/{run}_sorted_summits.bed", run=config["runs"].keys()),
        expand("analysis/peaks/{run}/{run}_treat_pileup.bw", run=config["runs"].keys()),
        expand("analysis/peaks/{run}/{run}_control_lambda.bw", run=config["runs"].keys()),
        #FASTQC_ALL
        expand("analysis/fastqc/{sample}_perSeqGC.txt", sample=config["samples"]),
        #CONSERVATION_ALL
        expand("analysis/peaks/{run}/{run}_sorted_5k_summits.bed", run=config["runs"].keys()),
        expand("analysis/conserv/{run}/{run}_conserv.txt", run=config["runs"].keys()),
        #CEAS_ALL
        expand("analysis/ceas/{run}/{run}_summary.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_stats.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_stats.txt", run=config["runs"].keys()),
        #FRIPS_ALL
        expand("analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam", sample=config["samples"]),
        expand("analysis/peaks/{run}/{run}_4M_peaks.narrowPeak", run=config["runs"].keys()),
        expand("analysis/frips/{run}/{run}_frip.txt",run=config["runs"].keys()),
        #MOTIF_ALL
        expand("analysis/motif/{run}/results/mdseqpos_out.html", run=config["runs"].keys()),
        #CONTAMINATION_ALL- need to figure how to handle this w/o doubling the code
        #expand("analysis/contam/{sample}/{sample}.{panel}.sai", sample=config['samples'].keys(), panel=_contaminationNames),
        #expand("analysis/contam/{sample}/{sample}.{panel}.bam", sample=config['samples'].keys(), panel=_contaminationNames),
        #expand("analysis/contam/{sample}/{sample}.{panel}.txt", sample=config['samples'].keys(), panel=_contaminationNames),
        #expand("analysis/contam/{sample}/{sample}_contamination.txt", sample=config['samples'].keys()),
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
