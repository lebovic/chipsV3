#MODULE: chilin - a chilin adapter of chips
# _logfile="analysis/logs/chilin.log"

import os

def chilin_targets(wildcards):
    ls = []
    
    for run in config["runs"]:
        ls.append("analysis/chilin/%s/%s_peaks.xls" % (run, run))
        if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
            ls.append("analysis/chilin/%s/%s_sorted_peaks.broadPeak.bed" % (run, run))
        else:
            ls.append("analysis/chilin/%s/%s_sorted_peaks.narrowPeak.bed" % (run, run))
            ls.append("analysis/chilin/%s/%s_sorted_summits.bed" % (run, run))
        ls.append("analysis/chilin/%s/%s_peaks.bed" % (run, run))
        ls.append("analysis/chilin/%s/%s_treat.bw" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_conserv_img.png" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_conserv.txt" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_gene_score_5fold.txt" % (run, run))
        ls.append("analysis/chilin/%s/attic/%s_gene_score.txt" % (run, run))
        if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] == False:
            if ("motif" in config) and config["motif"] == "mdseqpos":
                ls.append("analysis/chilin/%s/attic/%s_seqpos/" % (run, run))
        ls.append("analysis/chilin/%s/attic/json/" % run)
        #handle samples:
        run_samples = config['runs'][run]
        for sample in run_samples:
            if sample:
                ls.append("analysis/chilin/%s/attic/%s_100k_fastqc/" % (run, sample))
    return ls

def getJsonInput(wildcards):
    ls = []
    for run in config["runs"]:
        ls.append("analysis/json/%s/%s_conserv.json" % (run, run))
        if config["DHS"]:
            ls.append("analysis/json/%s/%s_dhs.json" % (run, run))        
        ls.append("analysis/json/%s/%s_frip.json" % (run, run))
        ls.append("analysis/json/%s/%s_macs2.json" % (run, run))
        if config['runs'][run][2]:
            ls.append("analysis/json/%s/%s_macs2_rep.json" % (run, run))
        ls.append("analysis/json/%s/%s_meta.json" % (run, run))
        ls.append("analysis/json/%s/%s_fastqc.json" % (run, run))
        ls.append("analysis/json/%s/%s_map.json" % (run, run))
        ls.append("analysis/json/%s/%s_enrich_meta.json" % (run, run))
        ls.append("analysis/json/%s/%s_pbc.json"% (run, run))
        ls.append("analysis/json/%s/%s_frag.json" % (run, run))
        if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] == False:
            if ("motif" in config) and config["motif"] == "mdseqpos":
                ls.append("analysis/json/%s/%s_seqpos.json" % (run, run))
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
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_sorted_peaks.narrowPeak.bed"
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getBroadPeakBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "analysis/chilin/{run}/{run}_sorted_peaks.broadPeak.bed"
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
       lambda wildcards: "analysis/motif/%s/results/mdseqpos_index.html" % chilin_getRunAndRep(wildcards)
   output:
       directory("analysis/chilin/{run}/attic/{run}_seqpos/")
   params:
       abspath = lambda wildcards: os.path.abspath(str("analysis/motif/%s" % chilin_getRunAndRep(wildcards)))
   shell:
       "ln -s {params.abspath}/* {output}"

rule getFastqc:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/"
    output:
        directory("analysis/chilin/{run}/attic/{sample}_100k_fastqc/")
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath}/* {output}"

rule getJson:
   input:
       getJsonInput
   output:
       directory("analysis/chilin/{run}/attic/json/")
   params:
       abspath = lambda wildcards, input: str(os.path.abspath(os.path.dirname(input[0])))
   shell:
       "ln -s {params.abspath}/* {output}"

#rule md5check:
#    input:
#        "analysis/chilin/{run}/"
#    output:
#        "analysis/chilin/{run}/{run}.md5"
#    shell:
#        "cidc_chips/modules/scripts/md5check.py -d {input} -I {wildcards.run}"
