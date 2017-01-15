#MODULE: qdnaseq- perform qdnaseq CNV analysis on sample bams
_logfile="analysis/logs/qdnaseq.log"

rule qdnaseq_all:
    input:
        "analysis/qdnaseq/qdnaseq.bed",
        "analysis/qdnaseq/qdnaseq.igv",
        "analysis/qdnaseq/qdnaseq.txt",
        "analysis/qdnaseq/qdnaseq.pdf",
        "analysis/qdnaseq/qdnaseq_segmented.igv",
        "analysis/qdnaseq/qdnaseq_calls.igv",
        #expand(".tmp/{sample}_unique.sorted.bam", sample=config["samples"]),

rule qdnaseq_linkFiles:
    """Link files into one single (temporary) directory"""
    input:
        "analysis/align/{sample}/{sample}_unique.sorted.bam"
    output:
        temp(".tmp/{sample}_unique.sorted.bam")
    message: "QDNAseq: linking files"
    log: _logfile
    shell:
        "ln -s ../{input} {output} 2>>{log}"

rule qdnaseq:
    """performs qdnaseq analysis on ALL of the samples (made by linkFiles)"""
    input:
        expand(".tmp/{sample}_unique.sorted.bam", sample=config["samples"]),
    output:
        "analysis/qdnaseq/qdnaseq.bed",
        "analysis/qdnaseq/qdnaseq.igv",
        "analysis/qdnaseq/qdnaseq.txt",
        "analysis/qdnaseq/qdnaseq.pdf",
        "analysis/qdnaseq/qdnaseq_segmented.igv",
        "analysis/qdnaseq/qdnaseq_calls.igv",
    message: "QDNAseq: performing cnv analysis"
    log: _logfile
    params:
        name="qdnaseq",
        qbin="chips/static/qdnaseq/qdnaseq_hg19_50.bin",
        out="analysis/qdnaseq/"
    shell:
        "R CMD BATCH --vanilla '--args {params.name} .tmp {params.qbin} {params.out}' chips/modules/scripts/qdnaseq.R {log}"

