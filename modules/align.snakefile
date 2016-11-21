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
        expand("analysis/align/{sample}/{sample}.sorted.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.sorted.bam", sample=config["samples"]),
        "analysis/align/mapping.csv",

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

rule uniquely_mapped_reads:
    """Get the uniquely mapped reads"""
    input:
        "analysis/align/{sample}/{sample}.bam"
    output:
        "analysis/align/{sample}/{sample}_unique.bam"
    message: "ALIGN: Filtering for uniquely mapped reads"
    log: _logfile
    shell:
        #NOTE: this is the generally accepted way of doing this as multiply 
        #mapped reads have a Quality score of 0
        "samtools view -bq 1 {input} > {output}"

rule map_stats:
    """Get the mapping stats for each aligment run"""
    input:
        bam="analysis/align/{sample}/{sample}.bam",
        uniq_bam="analysis/align/{sample}/{sample}_unique.bam"
    output:
        #temp("analysis/align/{sample}/{sample}_mapping.txt")
        "analysis/align/{sample}/{sample}_mapping.txt"
    message: "ALIGN: get mapping stats for each bam"
    log: _logfile
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "samtools flagstat {input.bam} > {output} 2>>{log}"
        " && samtools view -c {input.uniq_bam} >> {output} 2>> {log}"

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

rule sortBams:
    """General sort rule--take a bam {filename}.bam and 
    output {filename}.sorted.bam"""
    input:
        "analysis/align/{sample}/{filename}.bam"
    output:
        "analysis/align/{sample}/{filename}.sorted.bam"
    message: "ALIGN: sort bam file"
    log: _logfile
    threads: _bwa_threads
    shell:
        "samtools sort {input} -o {output} --threads {threads} 2>>{log}"
