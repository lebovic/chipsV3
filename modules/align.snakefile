#MODULE: Align fastq files to genome
#PARAMETERS:
_logfile="analysis/logs/align.log"
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"
_bwa_threads=4

rule align_all:
    input:
        expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"]),
        "analysis/align/mapping.csv"

def getFastq(wildcards):
    return config["samples"][wildcards.sample][int(wildcards.mate)]

def getFastq2(wildcards):
    return config["samples"][wildcards.sample]

def getMates(wildcards):
    s = wildcards.sample
    files = config["samples"][s]
    return ["analysis/align/%s/%s_%s.sai" % (s,s,m) for m in range(len(files))]

def getRunType(wildcards):
    s = wildcards.sample
    return "sampe" if len(config["samples"][s]) > 1 else "samse"


rule bwa_aln:
    input:
        getFastq
    output:
        sai="analysis/align/{sample}/{sample}_{mate}.sai"
    params:
        index=config['bwa_index'],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
    threads: _bwa_threads
    message: "ALIGN: Running BWA alignment"
    log: _logfile
    shell:
        "bwa aln -q {params.bwa_q} -l {params.bwa_l} -k {params.bwa_k} -t {threads} {params.index} {input} > {output.sai} 2>>{log}"

rule bwa_convert:
    input:
        sai=getMates,
        fastq=getFastq2
    output:
        "analysis/align/{sample}/{sample}.bam"
    params:
        run_type= getRunType,
        index=config['bwa_index'],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -"
    threads: _bwa_threads
    message: "ALIGN: Converting BWA alignment to BAM"
    log: _logfile
    shell:
        "bwa {params.run_type} {params.index} {input.sai} {input.fastq} | samtools {params.hack} > {output} 2>>{log}"

rule map_stats:
    """Get the mapping stats for each aligment run"""
    input:
        "analysis/align/{sample}/{sample}.bam"
    output:
        temp("analysis/align/{sample}/{sample}_mapping.txt")
        #"analysis/align/{sample}/{sample}_mapping.txt"
    message: "ALIGN: get mapping stats for each bam"
    log: _logfile
    shell:
        "samtools flagstat {input} > {output} 2>>{log}"

rule collect_map_stats:
    """Collect and parse out the mapping stats for the ALL of the samples"""
    input:
        expand("analysis/align/{sample}/{sample}_mapping.txt", sample=config["samples"])
    output:
        "analysis/align/mapping.csv"
    message: "ALIGN: collect and parse ALL mapping stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/align_getMapStats.py -f {files} > {output} 2>>{log}")
