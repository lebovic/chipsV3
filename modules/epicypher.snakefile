#MODULE: epicypher- quantify epicypher spike-ins
_logfile="analysis/logs/epicypher.log"
_threads=8

def epicypher_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/epicypher/%s/%s.epicypher.quant.txt" % (sample,sample))
    return ls

def getUnmappedReads(wildcards):
    sample = wildcards.sample
    ls = ["analysis/align/%s/%s.unmapped.fq.gz" % (sample, sample)]
    if len(config["samples"][wildcards.sample]) == 2:
        ls.append("analysis/align/%s/%s.unmapped.fq2.gz" % (sample, sample))
    return ls

rule epicypher_all:
    input:
        epicypher_targets

rule align_to_epicypher:
    """Align unmapped reads to epicypher assembly"""
    input:
        getUnmappedReads
    params:
        epicypher_index="chips/static/epicypher_bwa_index/epicypher.fa"
    output:
        "analysis/epicypher/{sample}/{sample}.epicypher.bam"
    threads: _threads
    message: "EPICYPHER: aligning unmapped reads to epicypher assembly"
    log: _logfile
    shell:
        "bwa mem -t {threads} {params.epicypher_index} {input} | samtools view -Sb - > {output}"

rule sort_epicypher:
    input:
        "analysis/epicypher/{sample}/{sample}.epicypher.bam"
    output:
        "analysis/epicypher/{sample}/{sample}.epicypher.sorted.bam"
    message: "EPICYPHER: sorting epicypher bam file"
    log: _logfile
    threads: _threads
    shell:
        "sambamba sort {input} -o {output} -t {threads} 2>>{log}"

rule uniquely_mapped:
    """Get uniquely mapped reads from epicypher.bam"""
    input:
        "analysis/epicypher/{sample}/{sample}.epicypher.sorted.bam"
    output:
        "analysis/epicypher/{sample}/{sample}.epicypher.sorted.unique.bam"
    message: "EPICYPHER: get uniquely mapped reads"
    threads: _threads
    log: _logfile
    shell:
        "samtools view -bq 1 -@ {threads} {input} > {output}"

rule index_epicypher:
    input:
        "analysis/epicypher/{sample}/{sample}.epicypher.sorted.unique.bam"
    output:
        "analysis/epicypher/{sample}/{sample}.epicypher.sorted.unique.bam.bai"
    message: "EPICYPHER: indexing epicypher bam file"
    log:_logfile
    shell:
        "sambamba index {input} {output}"

rule epicypher_quantify:
    """Quantify the reads for each spike-in mark"""
    input:
        bam="analysis/epicypher/{sample}/{sample}.epicypher.sorted.unique.bam",
        bai="analysis/epicypher/{sample}/{sample}.epicypher.sorted.unique.bam.bai"
    output:
        "analysis/epicypher/{sample}/{sample}.epicypher.quant.txt"
    message: "EPICYPHER: quantifying epicypher marks"
    threads: _threads
    log: _logfile
    shell:
        "chips/modules/scripts/epicypher_quant.py -b {input.bam} > {output}"
