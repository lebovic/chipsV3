#MODULE: motif module--uses MDSeqPos.py to perform motif analysis and generate 
#motif table
_logfile="analysis/logs/motif.log"

rule motif_all:
    input:
        expand("analysis/motif/{run}/results/mdseqpos_out.html", run=config["runs"].keys()),
        #expand("analysis/conserv/{run}/{run}_conserv.txt", run=config["runs"].keys()),

rule motif:
    """call MDSeqpos on top 5k summits"""
    input:
        "analysis/peaks/{run}/{run}_sorted_5k_summits.bed"
    output:
        results="analysis/motif/{run}/results",
        html="analysis/motif/{run}/results/mdseqpos_out.html",
    params:
        path=config['motif_path'],
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        run = lambda wildcards: wildcards.run,
    message: "MOTIF: calling MDSeqPos on top 5k summits"
    log: _logfile
    shell:
        "{params.pypath} {config[mdseqpos_path]} {input} {params.path} -m cistrome.xml -d -O analysis/motif/{params.run}/results 1>>{log}"
