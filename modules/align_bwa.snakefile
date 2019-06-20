#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
# _logfile="analysis/logs/align.log"
_bwa_threads=8
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"

def getFastq(wildcards):
    return config["samples"][wildcards.sample]


def getAlnFastq(wildcards):
    return config["samples"][wildcards.sample][int(wildcards.mate)]


def getMates(wildcards):
    s = wildcards.sample
    files = config["samples"][s]
    return ["analysis/align/%s/%s_%s_aln.sai" % (s,s,m) for m in range(len(files))]


def getRunType(wildcards):
    s = wildcards.sample
    return "sampe" if len(config["samples"][s]) > 1 else "samse"


checkpoint read_length:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        # direc=directory("analysis/fastqc/{sample}/"),
        length=temp("analysis/fastqc/{sample}/{sample}_read_length.txt")
    shell:
        # "mkdir analysis/fastqc/{wildcards.sample}/;"
        "python cidc_chips/modules/scripts/align_get_read_length.py -f {input} -o {output.length}"


rule bwa_mem:
    input:
        getFastq
    output:
        temp("analysis/align/{sample}/{sample}_mem.bam")
    params:
        index=config['bwa_index'],
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else ""
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment"
    log: "analysis/logs/align/{sample}.log"
    conda: "../envs/align/align_bwa.yaml"
    shell:
        "{params.sentieon} bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"


rule bwa_aln:
    input:
        getAlnFastq
    output:
        sai=temp("analysis/align/{sample}/{sample}_{mate}_aln.sai")
    params:
        index=config['bwa_index'],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else ""
    threads: _bwa_threads
    message: "ALIGN: Running BWA alignment"
    # log: "analysis/logs/align/{sample}.log"
    shell:
        "{params.sentieon} bwa aln -q {params.bwa_q} -l {params.bwa_l} -k {params.bwa_k} -t {threads} {params.index} {input} > {output.sai}" # 2>>{log}"

rule bwa_convert:
    input:
        sai=getMates,
        fastq=getFastq
    output:
        temp("analysis/align/{sample}/{sample}_aln.bam")
    params:
        run_type= getRunType,
        index=config['bwa_index'],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -",
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else ""
    threads: _bwa_threads
    message: "ALIGN: Converting BWA alignment to BAM"
    # log: "analysis/logs/align/{sample}.log"
    shell:
        "{params.sentieon} bwa {params.run_type} {params.index} {input.sai} {input.fastq} | samtools {params.hack} > {output}" # 2>>{log}"


def aggregate_align_input(wildcards):
    # decision based on content of output file
    with open(checkpoints.read_length.get(sample=wildcards.sample).output[0]) as f:
        if int(f.read().strip()) >= 40:
            return "analysis/align/{sample}/{sample}_mem.bam"
        else:
            return "analysis/align/{sample}/{sample}_aln.bam"


rule aggregate_align:
    input:
        aggregate_align_input
    output:
        temp("analysis/align/{sample}/{sample}.bam")
    shell:
        "mv {input} {output}"
