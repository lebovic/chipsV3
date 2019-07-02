#MODULE: chilin - a chilin adapter of chips
# _logfile="analysis/logs/chilin.log"
import os

if "Chilin_path" in config and config["Chilin_path"]:
    chilinpath = config["Chilin_path"]
else:
    chilinpath = "analysis/chilin"

def chilin_targets(wildcards):
    ls = []

    for run in config["runs"]:
        ls.append("%s/%s/%s_peaks.xls" % (chilinpath, run, run))
        if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
            ls.append("%s/%s/%s_sorted_peaks.broadPeak.bed" % (chilinpath, run, run))
        else:
            ls.append("%s/%s/%s_sorted_peaks.narrowPeak.bed" % (chilinpath, run, run))
            ls.append("%s/%s/%s_sorted_summits.bed" % (chilinpath, run, run))
        ls.append("%s/%s/%s_peaks.bed" % (chilinpath, run, run))
        ls.append("%s/%s/%s_treat.bw" % (chilinpath, run, run))
        ls.append("%s/%s/attic/%s_conserv_img.png" % (chilinpath, run, run))
        ls.append("%s/%s/attic/%s_conserv.txt" % (chilinpath, run, run))
        ls.append("%s/%s/attic/%s_gene_score_5fold.txt" % (chilinpath, run, run))
        ls.append("%s/%s/attic/%s_gene_score.txt" % (chilinpath, run, run))
        if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] != True:
            if ("motif" in config) and config["motif"] == "mdseqpos":
                ls.append("%s/%s/attic/%s_seqpos/" % (chilinpath, run, run))
        ls.append("%s/%s/attic/json/" % (chilinpath, run))
        #handle samples:
        run_samples = config['runs'][run]
        for sample in run_samples:
            if sample:
                ls.append("%s/%s/attic/%s_100k_fastqc/" % (chilinpath, run, sample))
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
        "%s/{run}/{run}_peaks.xls" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getNarrowPeakBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/{run}_sorted_peaks.narrowPeak.bed" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getBroadPeakBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/{run}_sorted_peaks.broadPeak.bed" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getSortedSummitsBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_sorted_summits.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/{run}_sorted_summits.bed" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getPeaksBed:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.bed" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/{run}_peaks.bed" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getTreatBw:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_treat_pileup.bw" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/{run}_treat.bw" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservPng:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.png" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/attic/{run}_conserv_img.png" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getConservTxt:
    input:
        lambda wildcards: "analysis/conserv/%s/%s_conserv.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/attic/{run}_conserv.txt" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getReg5foldTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score_5fold.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/attic/{run}_gene_score_5fold.txt" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getRegTxt:
    input:
        lambda wildcards: "analysis/regulatory/%s/%s_gene_score.txt" % (chilin_getRunAndRep(wildcards), chilin_getRunAndRep(wildcards))
    output:
        "%s/{run}/attic/{run}_gene_score.txt" % chilinpath
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath} {output}"

rule getMotif:
   input:
       lambda wildcards: "analysis/motif/%s/results/mdseqpos_index.html" % chilin_getRunAndRep(wildcards)
   output:
       directory("%s/{run}/attic/{run}_seqpos/" % chilinpath )
   params:
       abspath = lambda wildcards: os.path.abspath(str("analysis/motif/%s" % chilin_getRunAndRep(wildcards)))
   shell:
       "ln -s {params.abspath}/* {output}"

rule getFastqc:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/"
    output:
        directory("%s/{run}/attic/{sample}_100k_fastqc/" % chilinpath)
    params:
        abspath = lambda wildcards, input: os.path.abspath(str(input))
    shell:
        "ln -s {params.abspath}/* {output}"

rule getJson:
   input:
       getJsonInput
   output:
       directory("%s/{run}/attic/json/" % chilinpath)
   params:
       abspath = lambda wildcards, input: str(os.path.abspath(os.path.dirname(input[0])))
   shell:
       "ln -s {params.abspath}/* {output}"

#rule md5check:
#    input:
#        "%s/{run}/"
#    output:
#        "%s/{run}/{run}.md5 % chilinpath"
#    shell:
#        "cidc_chips/modules/scripts/md5check.py -d {input} -I {wildcards.run}"
