#MODULE: motif module--uses MDSeqPos.py to perform motif analysis and generate 
#motif table
_logfile="analysis/logs/motif.log"

def motif_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        ls.append("analysis/motif/%s/results/mdseqpos_index.html" % run)
        ls.append("analysis/motif/%s/results/motif_list.json" % run)
    return ls

rule motif_all:
    input:
        motif_targets

rule motif:
    """call MDSeqpos on top 5k summits"""
    input:
        "analysis/peaks/{run}/{run}_sorted_5k_summits.bed"
    output:
        results="analysis/motif/{run}/results",
        #html="analysis/motif/{run}/results/mdseqpos_out.html",
        html="analysis/motif/{run}/results/mdseqpos_index.html",
        json="analysis/motif/{run}/results/motif_list.json",
    params:
        path=config['motif_path'],
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        run = lambda wildcards: wildcards.run,
    message: "MOTIF: calling MDSeqPos on top 5k summits"
    log: _logfile
    shell:
        "{params.pypath} {config[mdseqpos_path]} {input} {params.path} -m cistrome.xml -d -O analysis/motif/{params.run}/results 1>>{log}"
