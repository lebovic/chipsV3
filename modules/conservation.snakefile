#MODULE: conservation- module to create conservation plots
_logfile="analysis/logs/conservation.log"

def conservation_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        ls.append("analysis/peaks/%s/%s_sorted_5k_summits.bed" % (run,run))
        ls.append("analysis/conserv/%s/%s_conserv.R" % (run,run))
        ls.append("analysis/conserv/%s/%s_conserv.pdf" % (run,run))
    ls.append("analysis/conserv/conservationPlot.png")
    return ls

rule conservation_all:
    input:
        conservation_targets

rule top5k_peaks:
    """take the top 5000 peaks, sorted by score"""
    input:
        "analysis/peaks/{run}/{run}_sorted_summits.bed"
    output:
        "analysis/peaks/{run}/{run}_sorted_5k_summits.bed"
    params:
        lines = 5000
    message: "CONSERVATION: top5k_peaks"
    log: _logfile
    shell:
        "head -n {params.lines} {input} > {output} 2>>{log}"

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
        run = lambda wildcards: wildcards.run,
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        #name= wildcards.run
    message: "CONSERVATION: calling conservation script"
    log: _logfile
    shell:
        "{params.pypath} {config[python2]} chips/modules/scripts/conservation_plot.py -t Conservation_at_summits -d {params.db} -o analysis/conserv/{params.run}/{params.run}_conserv -l Peak_summits {input} -w {params.width} > {output.score} 2>>{log}"

rule conservation_plot:
    """generate (overall) conservation plot for ALL runs"""
    input:
        expand("analysis/conserv/{run}/{run}_conserv.R", run=config['runs'])
    output:
        png="analysis/conserv/conservationPlot.png",
    params:
        rout = "analysis/conserv/conservationPlot.R",
        template = "chips/static/conserv_plotConserv.R.txt"
    message: "CONSERVATION: generating OVERALL conservation plot"
    log: _logfile
    run:
        files = " -r ".join(input)
        #Generate and execute analysis/conserv/conservationPlot.R
        shell("chips/modules/scripts/conserv_plotConserv.py -r {files} -t {params.template} -p {output.png} -o {params.rout} && Rscript {params.rout} 2>>{log}")

