
#TODO: handle control
def getTreats(wildcards):
    r = config['runs'][wildcards.run]
    #convert SAMPLE names to BAMS
    tmp = ["analysis/align/%s/%s.bam" % (s,s) for s in r if s]
    #print(tmp)
    return tmp[:2]

def getConts(wildcards):
    r = config['runs'][wildcards.run]
    return r[2:4]

rule peaks_all:
    input:
        expand("analysis/peaks/{run}/{run}_sorted_summits.bed", run=config["runs"].keys()),

rule macs2_callpeaks:
    input:
        treat=getTreats
    output:
        "analysis/peaks/{run}/{run}_peaks.narrowPeak",
        #"analysis/peaks/{run}/{run}_peaks.bed",
        "analysis/peaks/{run}/{run}_peaks.xls",
        "analysis/peaks/{run}/{run}_summits.bed",
        "analysis/peaks/{run}/{run}_treat_pileup.bdg",
        "analysis/peaks/{run}/{run}_control_lambda.bdg",
    params:
        fdr="0.01",
        extsize="146",
        species="hs",
        outdir="analysis/peaks/{run}/",
        name="{run}"
    shell:
        "macs2 callpeak --SPMR -B -q {params.fdr} --keep-dup 1 -g {params.species} --extsize {params.extsize} --nomodel -t {input.treat} --outdir {params.outdir} -n {params.name}"


rule peakToBed:
    input:
        "analysis/peaks/{run}/{run}_peaks.narrowPeak"
    output:
        "analysis/peaks/{run}/{run}_peaks.bed"
    message: "PEAKS: Converting peak file to bed file"
    shell:
        "cut -f1,2,3,4,9 {input} > {output}"

rule sortSummits:
    input:
        "analysis/peaks/{run}/{run}_summits.bed"
    output:
        "analysis/peaks/{run}/{run}_sorted_summits.bed"
    message: "PEAKS: sorting the summits bed by score"
    shell:
        "sort -r -n -k 5 {input} > {output}"
