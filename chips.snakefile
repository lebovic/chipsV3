#!/usr/bin/env python

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
#-----------------------------------------

rule target:
    #input: expand("analysis/{sample}/{sample}.sai", sample=config["samples"].keys())
    #input: expand("analysis/{sample}/{sample}.sai", sample=config["samples"])
    input: expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"])
    message: "Compiling all output"

include: "./modules/align.snakefile"         # rules specific to BWA

    
#TODO:
#1. handle PE--which is bwa sampe [index] sai1 sai2 fastq1 fastq2
#2. peak calling
