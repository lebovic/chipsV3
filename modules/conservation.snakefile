#MODULE: conservation- module to create conservation plots

rule conservation_all:
    input:
        expand("analysis/peaks/{run}/{run}_sorted_5k_summits.bed", run=config["runs"].keys()),
        expand("analysis/conserv/{run}/{run}_conserv.txt", run=config["runs"].keys()),

rule top5k_peaks:
    """take the top 5000 peaks, sorted by score"""
    input:
        "analysis/peaks/{run}/{run}_sorted_summits.bed"
    output:
        "analysis/peaks/{run}/{run}_sorted_5k_summits.bed"
    params:
        lines = 5000
    shell:
        "head -n {params.lines} {input} > {output}"

rule conservation:
    """generate conservation plots"""
    input:
        "analysis/peaks/{run}/{run}_sorted_summits.bed"
    output:
        pdf="analysis/conserv/{run}/{run}_conserv.pdf",
        r="analysis/conserv/{run}/{run}_conserv.R",
        score="analysis/conserv/{run}/{run}_conserv.txt",
    params:
        db=config['conservation'],
        width=4000,
        run = lambda wildcards: wildcards.run
        #name= wildcards.run
    log: "log.txt"
    shell:
        "conserv/conservation_plot.py -t Conservation_at_summits -d {params.db} -o analysis/conserv/{params.run}/{params.run}_conserv -l Peak_summits {input} -w {params.width} > {output.score}"
