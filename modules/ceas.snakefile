#MODULE: CEAS- annotating where the peaks fall (promoter, exon, intron, interg)
_logfile="analysis/logs/ceas.log"

rule ceas_all:
    input:
        expand("analysis/ceas/{run}/{run}_summary.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_DHS_stats.txt", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_peaks.bed", run=config["runs"].keys()),
        expand("analysis/ceas/{run}/{run}_velcro_stats.txt", run=config["runs"].keys()),

rule ceas:
    """Annotate peak regions"""
    input:
        "analysis/peaks/{run}/{run}_peaks.bed"
    output:
        promoter="analysis/ceas/{run}/{run}_peaks_promoter.bed",
        exon="analysis/ceas/{run}/{run}_peaks_exon.bed",
        intron="analysis/ceas/{run}/{run}_peaks_intron.bed",
        intergenic="analysis/ceas/{run}/{run}_peaks_intergenic.bed",
        summary="analysis/ceas/{run}/{run}_summary.txt",
    message: "CEAS: annotating peak regions"
    log: _logfile
    params:
        db=config['geneTable'],
        path="analysis/ceas/{run}/"
    shell:
        "chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} > {output.summary} 2>>{log}"


rule takeTop5k:
    """Take the top 5000 sites"""
    input:
        "analysis/peaks/{run}/{run}_sorted_peaks.bed"
    params:
        n=5000
    message: "DHS: Take top sites"
    log: _logfile
    output:
        temp('analysis/ceas/{run}/{run}_sorted_5k_peaks.bed')
    shell:
        "head -n {params.n} {input} > {output} 2>>{log}"

#------------------------------------------------------------------------------
rule DHS_intersectDHS:
    """Intersect PEAKS with DHS regions"""
    input:
        'analysis/ceas/{run}/{run}_sorted_5k_peaks.bed'
    params:
        dhs=config['DHS']
    message: "DHS: intersect PEAKS with DHS regions"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_DHS_peaks.bed'
    shell:
        "intersectBed -wa -u -a {input} -b {params.dhs} > {output} 2>>{log}"

rule DHS_stat:
    """collect DHS stats"""
    input:
        n='analysis/ceas/{run}/{run}_sorted_5k_peaks.bed',
        dhs='analysis/ceas/{run}/{run}_DHS_peaks.bed'
    message: "DHS: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_DHS_stats.txt'
    shell:
        "wc -l {input.n} {input.dhs} > {output} 2>>{log}"
#------------------------------------------------------------------------------
rule VELCRO_intersectVelcro:
    """Intersect PEAKS with velcro regions"""
    input:
        'analysis/ceas/{run}/{run}_sorted_5k_peaks.bed'
    params:
        velcro=config['velcro_regions']
    message: "VELCRO: intersect PEAKS with velcro regions"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_velcro_peaks.bed'
    shell:
        "intersectBed -wa -u -a {input} -b {params.velcro} > {output} 2>>{log}"
    
rule VELCRO_stat:
    """collect VELCRO stats"""
    input:
        n='analysis/ceas/{run}/{run}_sorted_5k_peaks.bed',
        velcro='analysis/ceas/{run}/{run}_velcro_peaks.bed'
    message: "VELCRO: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_velcro_stats.txt'
    shell:
        "wc -l {input.n} {input.velcro} > {output} 2>>{log}"
