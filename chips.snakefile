#!/usr/bin/env python

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
#-----------------------------------------

rule target:
    input: expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"])
    message: "Compiling all output"

include: "./modules/align.snakefile"         # rules specific to BWA
#include: "./modules/peaks.snakefile"         # peak calling rules
