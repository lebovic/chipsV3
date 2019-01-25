#MODULE: chilin - a chilin adapter of chips
_logfile="analysis/logs/chilin.log"

import os

def chilin_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/chilin/%s/%s_peaks.xls" % (sample, sample))
        ls.append("analysis/chilin/%s/%s_sorted_peaks.narrowPeak.bed" % (sample, sample))
        ls.append("analysis/chilin/%s/%s_sorted_summits.bed" % (sample, sample))
        ls.append("analysis/chilin/%s/%s_peaks.bed" % (sample, sample))
        ls.append("analysis/chilin/%s/%s_treat.bw" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_conserv_img.png" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_conserv.txt" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_gene_score_5fold.txt" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_gene_score.txt" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_motif/" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/%s_100k_fastqc/" % (sample, sample))
        ls.append("analysis/chilin/%s/attic/json/" % sample)
        ls.append("analysis/chilin/%s/%s.md5" % (sample, sample))
    return ls

def getRunAndRep(wildcards):
    item = 0
    replicate = []
    for i in range(2,len(r)+2,2):
        a = int(i/2)
        replicate.append("rep%s" % str(a))
        replicate.append("rep%s" % str(a))
        item += 2
    # print(replicate)
    sample = wildcards.sample
    # print(sample)
    key_list=[]
    value_list=[]
    for key,value in config["runs"].items():
        key_list.append(key)
        value_list.append(value)
    for i in value_list:
        if sample in i:
            run = key_list[value_list.index(i)]
            rep = replicate[i.index(sample)]
            break
        else:
            continue
    return "%s.%s" % (run, rep)


rule chilin_all:
    input:
        chilin_targets

rule getPeaksXls:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.xls" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/{sample}_peaks.xls"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getNarrowPeakBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.narrowPeak.bed" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/{sample}_sorted_peaks.narrowPeak.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getSortedSummitsBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_summits.bed" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/{sample}_sorted_summits.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getPeaksBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.bed" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/{sample}_peaks.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getTreatBw:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_treat_pileup.bw" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/{sample}_treat.bw"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservPng:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.png" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/attic/{sample}_conserv_img.png"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservTxt:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.txt" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/attic/{sample}_conserv.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getReg5foldTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score_5fold.txt" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/attic/{sample}_gene_score_5fold.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getRegTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score.txt" % (getRunAndRep(wildcards), getRunAndRep(wildcards))
    output:
        "analysis/chilin/{sample}/attic/{sample}_gene_score.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getMotif:
    input:
        lambda wildcards: "analysis/motif/%s/" % getRunAndRep(wildcards)
    output:
        "analysis/chilin/{sample}/attic/{sample}_motif/"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getFastqc:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/"
    output:
        "analysis/chilin/{sample}/attic/{sample}_100k_fastqc/"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getJson:
    input:
        "analysis/json/{sample}/"
    output:
        "analysis/chilin/{sample}/attic/json/"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule md5check:
    input:
        "analysis/chilin/{sample}/"
    output:
        "analysis/chilin/{sample}/{sample}.md5"
    shell:
        "cidc_chips/modules/scripts/md5check.py -d {input} -I {wildcards.sample}"



