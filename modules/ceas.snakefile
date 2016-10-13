#MODULE: CEAS- annotating where the peaks fall (promoter, exon, intron, interg)

rule ceas_all:
    input:
        expand("analysis/ceas/{run}/{run}_summary.txt", run=config["runs"].keys()),

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
    params:
        db=config['geneTable'],
        path="analysis/ceas/{run}/"
    shell:
        "chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} > {output.summary}"
