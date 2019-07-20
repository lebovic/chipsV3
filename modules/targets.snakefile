#MODULE: targets module--use BETA to calculate the Regular Potential

#PARAMETERS
# _logfile="analysis/logs/targets.log"

def targets_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/targets/%s/%s_gene_score_5fold.txt" % (runRep,runRep))
            ls.append("analysis/targets/%s/%s_gene_score.txt" % (runRep,runRep))
    return ls

rule targets_all:
    input:
        targets_targets

rule get_5Fold_Peaks:
    input:
        # "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.xls"
        "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.bed"
    output:
        "analysis/targets/{run}.{rep}/{run}.{rep}_5foldpeaks.bed"
    message: "REGULATORY: get 5 fold peaks"
    log:"analysis/logs/targets/{run}.{rep}.log"
    shell:
        # """awk '($1 != "chr" && $1 !="#" && $8 >= 5)' {input} | awk '{{OFS="\\t"; print $1,$2,$3,$10,$9}}' > {output}"""
        "awk '($5 >= 5)' {input} > {output}"

rule get_5Fold_Peaks_RP_Score:
    input:
        "analysis/targets/{run}.{rep}/{run}.{rep}_5foldpeaks.bed"
    output:
        "analysis/targets/{run}.{rep}/{run}.{rep}_gene_score_5fold.txt"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        genome=config['geneTable']
    message: "REGULATORY: get RP score of 5 fold peaks"
    log:"analysis/logs/targets/{run}.{rep}.log"
    shell:
        "{params.pypath} cidc_chips/modules/scripts/targets_getRP.py -t {input} -g {params.genome} -n {output} -d 100000"

rule get_top_peaks:
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        "analysis/targets/{run}.{rep}/{run}.{rep}_peaks_top_reg.bed"
    params:
        peaks = 10000
    log:"analysis/logs/targets/{run}.{rep}.log"
    message: "REGULATORY: get top summits for regpotential"
    shell:
        "head -n {params.peaks} {input} > {output}"

rule get_top_Peaks_RP_Score:
    input:
        "analysis/targets/{run}.{rep}/{run}.{rep}_peaks_top_reg.bed"
    output:
        "analysis/targets/{run}.{rep}/{run}.{rep}_gene_score.txt"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        genome=config['geneTable']
    message: "REGULATORY: get RP score of top peaks"
    log:"analysis/logs/targets/{run}.{rep}.log"
    shell:
        "{params.pypath} cidc_chips/modules/scripts/targets_getRP.py -t {input} -g {params.genome} -n {output} -d 100000"

