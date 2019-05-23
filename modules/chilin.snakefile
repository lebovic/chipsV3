#MODULE: chilin - a chilin adapter of chips
# _logfile="analysis/logs/chilin.log"

import os

def chilin_targets(wildcards):
    ls = []
    
    for run in config["runs"]:
        ls.append("analysis/chilin/%s/%s_peaks.xls" % (run, run))
        ls.append("analysis/chilin/%s/%s_sorted_peaks.narrowPeak.bed" % (run, run))
        ls.append("analysis/chilin/%s/%s_sorted_summits.bed" % (run, run))
        ls.append("analysis/chilin/%s/%s_peaks.bed" % (run, run))
        ls.append("analysis/chilin/%s/%s_treat.bw" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_conserv_img.png" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_conserv.txt" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_gene_score_5fold.txt" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_gene_score.txt" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_motif/" % (run, run))
        ls.append("analysis/chilin/%s/attic/json/" % run)
        #handle samples:
        run_samples = config['runs'][run]
        for sample in run_samples:
            if sample:
                ls.append("analysis/chilin/%s/attic/%s_100k_fastqc/" % (run, sample))

    return ls

def chilin_getRunAndRep(wildcards):
    #Xindong, this is a hack until we figure out a cleaner way to do this
    #PROBLEM: the chilin "run" does not have any idea about replicates
    #chips runs do.
    #HACK/Solution: we just eliminate the rep part because cistromedb
    #does not run replicates; we assume chips's {run}.{rep1} = chilin {run}

    run_rep = "%s.rep1" % wildcards.run
    return run_rep
    # item = 0
    # replicate = []
    # for i in range(2,len(r)+2,2):
    #     a = int(i/2)
    #     replicate.append("rep%s" % str(a))
    #     replicate.append("rep%s" % str(a))
    #     item += 2
    # # print(replicate)
    # sample = wildcards.sample
    # # print(sample)
    # key_list=[]
    # value_list=[]
    # for key,value in config["runs"].items():
    #     key_list.append(key)
    #     value_list.append(value)
    # for i in value_list:
    #     if sample in i:
    #         run = key_list[value_list.index(i)]
    #         rep = replicate[i.index(sample)]
    #         break
    #     else:
    #         continue
    # return "%s.%s" % (run, rep)

rule chilin_all:
    input:
        chilin_targets

rule getPeaksXls:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.xls" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_peaks.xls"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getNarrowPeakBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.narrowPeak.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_sorted_peaks.narrowPeak.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getSortedSummitsBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_summits.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_sorted_summits.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getPeaksBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_peaks.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getTreatBw:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_treat_pileup.bw" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_treat.bw"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservPng:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.png" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/attic/{run}_conserv_img.png"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservTxt:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/attic/{run}_conserv.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getReg5foldTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score_5fold.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/attic/{run}_gene_score_5fold.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getRegTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/attic/{run}_gene_score.txt"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getMotif:
   input:
       lambda wildcards: "analysis/motif/%s/" % chilin_getRunAndRep(wildcards)
   output:
       "analysis/chilin/{run}/attic/{run}_motif/"
   params:
       abspath = lambda wildcards, input: os.path.abspath(str(input))
   shell:
       "ln -s {params.abspath}/* {output}"

rule getFastqc:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/"
    output:
        "analysis/chilin/{run}/attic/{sample}_100k_fastqc/"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath}/* {output}"

rule getJson:
   input:
       "analysis/json/{run}/"
   output:
       "analysis/chilin/{run}/attic/json/"
   params:
       abspath = lambda wildcards, input: os.path.abspath(str(input))
   shell:
       "ln -s {params.abspath}/* {output}"

#rule md5check:
#    input:
#        "analysis/chilin/{run}/"
#    output:
#        "analysis/chilin/{run}/{run}.md5"
#    shell:
#        "cidc_chips/modules/scripts/md5check.py -d {input} -I {wildcards.run}"
