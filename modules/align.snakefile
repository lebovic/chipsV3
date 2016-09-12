#MODULE: Align fastq files to genome

#rule bwa_all:
#    input:
#DICTIONARY of whether samples are SE or PE
_mates = dict([(s, len(config["samples"][s])) for s in config["samples"]])
print(_mates)

def getFastq(wildcards):
    return config["samples"][wildcards.sample][int(wildcards.mate)]

def getFastq2(wildcards):
    return config["samples"][wildcards.sample]

def getMates(wildcards):
    s = wildcards.sample
    files = config["samples"][s]
    return ["analysis/align/%s/%s_%s.sai" % (s,s,m) for m in range(len(files))]

rule bwa_aln:
    input:
        getFastq
    output:
        sai="analysis/align/{sample}/{sample}_{mate}.sai"
    params:
        index=config['bwa_index'],
    threads: 4
    message: "ALIGN: Running BWA alignment"
    shell:
        "bwa aln -q 5 -l 32 -k 2 -t {threads} {params.index} {input} > {output.sai}"

rule bwa_convert:
    input:
        #sai="analysis/align/{sample}/{sample}_{mate}.sai",
        sai=getMates,
        fastq=getFastq2
    output:
        "analysis/align/{sample}/{sample}.bam"
    params:
        index=config['bwa_index'],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -"
    threads: 4
    message: "ALIGN: Converting BWA alignment to BAM"
    shell:
        "bwa samse {params.index} {input.sai} {input.fastq} | samtools {params.hack} > {output}"
