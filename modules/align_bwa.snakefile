#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
# _logfile=output_path + "/logs/align.log"
import subprocess

_bwa_threads=8
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"

_bwa_ram=7 #7GB

def getFastq(wildcards):
    return config["samples"][wildcards.sample]


def getAlnFastq(wildcards):
    return config["samples"][wildcards.sample][int(wildcards.mate)]

def getBam(wildcards):
    return config["samples"][wildcards.sample][0] #returns only the 1st elm

def getMates(wildcards):
    s = wildcards.sample
    files = config["samples"][s]
    return ["%s/align/%s/%s_%s_aln.sai" % (output_path,s,s,m) for m in range(len(files))]


def getRunType(wildcards):
    s = wildcards.sample
    return "sampe" if len(config["samples"][s]) > 1 else "samse"


checkpoint align_readsLength:
    input:
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        # direc=directory(output_path + "/fastqc/{sample}/"),
        length=temp(output_path + "/fastqc/{sample}/{sample}_read_length.txt")
    shell:
        # "mkdiroutput_path +  /fastqc/{wildcards.sample}/;"
        src_path + "/modules/scripts/align_getReadLength.py -f {input} -o {output.length}"


rule align_bwaMem:
    input:
        getFastq
    output:
        temp(output_path + "/align/{sample}/{sample}_mem.bam")
    params:
        index=config['bwa_index'],
        sentieon=config.get("sentieon", ""),
        read_group=lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample)
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment for {input}"
    log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/benchmark/{sample}_align_bwaMem.benchmark"
    conda: "../envs/align/align_bwa.yaml"
    resources:
        mem_mb = 1024*_bwa_ram,
    shell:
        "{params.sentieon} bwa mem -t {threads} -R \"{params.read_group}\" {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"


rule align_bwaAln:
    input:
        getAlnFastq
    output:
        sai=temp(output_path + "/align/{sample}/{sample}_{mate}_aln.sai")
    params:
        index=config['bwa_index'],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else ""
    threads: _bwa_threads
    message: "ALIGN: Running BWA aln for alignment for {input}"
    # log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/benchmark/{sample}_{mate}_align_bwaAln.benchmark"
    conda: "../envs/align/align_bwa.yaml"
    resources:
        mem_mb = 1024*_bwa_ram,
    shell:
        "{params.sentieon} bwa aln -q {params.bwa_q} -l {params.bwa_l} -k {params.bwa_k} -t {threads} {params.index} {input} > {output.sai}"

rule align_bwaConvert:
    input:
        sai=getMates,
        fastq=getFastq
    output:
        temp(output_path + "/align/{sample}/{sample}_aln.bam")
    params:
        run_type= getRunType,
        index=config['bwa_index'],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -",
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else "",
        read_group=lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample)
    threads: _bwa_threads
    message: "ALIGN: Converting BWA alignment to BAM"
    # log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/benchmark/{sample}_align_bwaConvert.benchmark"
    conda: "../envs/align/align_bwa.yaml"
    resources:
        mem_mb = 1024*_bwa_ram,
    shell:
        """{params.sentieon} bwa {params.run_type} -r \"{params.read_group}\" {params.index} {input.sai} {input.fastq} | samtools {params.hack} > {output}"""

rule align_from_bam:
    input:
        getBam
    output:
        bam="analysis/align/{sample}/{sample}_fromBam.bam",
    threads: _bwa_threads
    params:
        sentieon=config.get('sentieon',''),
        index=config['bwa_index'],
        #DON'T write a mini-program withi a program-
        #awk cmd to add sample names to RGs!!
        awk_cmd=lambda wildcards: "awk -v OFS=\'\\t\' \'{ split($2,a,\":\"); read_id=a[2]; $2=\"ID:%s.\" read_id; gsub(/SM:.+\\t/,\"SM:%s\\t\"); print $0}\'" % (wildcards.sample, wildcards.sample),
        #NEVER do it twice!- gawk cmd to inject sample name into each read!!!
        gawk_cmd=lambda wildcards: "gawk -v OFS=\'\\t\' \'{rg=match($0,/RG:Z:(\S+)/,a); read_id=a[1]; if (rg) {sub(/RG:Z:\S+/, \"RG:Z:%s.\" read_id, $0); print $0} else { print $0 }}\'" % wildcards.sample,
    benchmark: output_path + "/benchmark/align/{sample}/{sample}.align_from_bam.txt"
    conda: "../envs/align/align_bwa.yaml"
    resources:
        mem_mb = 1024*_bwa_ram,
    shell:
        """samtools view -H {input} | grep \"^@RG\" | {params.awk_cmd} > {wildcards.sample}.header && samtools collate --output-fmt SAM -@ {threads} -Of {input} | {params.gawk_cmd} | samtools view -@ {threads} -b - | samtools fastq -@ {threads} -t -s /dev/null -0 /dev/null - | ({params.sentieon} bwa mem -t {threads} -M -K 10000000 -p -C -H {wildcards.sample}.header {params.index} - || echo -n 'error' ) | samtools view -Sb - > {output}; rm {wildcards.sample}.header"""

def aggregate_align_input(wildcards):
    # decision based on 1. whether it's a bam or fastq
    # IF it's a fastq, then decide whether it's >= 40 bp (bwa mem) or bwa aln
    sample_first_file = config["samples"][wildcards.sample][0]
    ret = ""
    if sample_first_file.endswith('.bam'):
        ret= output_path + "/align/{sample}/{sample}_fromBam.bam"
    else: #it's a fastq
        with open(checkpoints.align_readsLength.get(sample=wildcards.sample).output[0]) as f:
            if int(f.read().strip()) >= 40:
                ret= output_path + "/align/{sample}/{sample}_mem.bam"
            else:
                ret= output_path + "/align/{sample}/{sample}_aln.bam"
    #print(ret)
    return ret

checkpoint align_aggregate:
    input:
        aggregate_align_input
    output:
        temp(output_path + "/align/{sample}/{sample}.bam")
    shell:
        "mv {input} {output}"


rule align_macsRunInfo:
    """Dump the current version of bwa into a text file for the report"""
    output:
        output_path + "/align/run_info.txt"
    message: "ALIGN/REPORT - collection bwa version info"
    conda: "../envs/align/align_bwa.yaml"
    shell:
        src_path + "/modules/scripts/align_parseBwaVersion.py -o {output}"
