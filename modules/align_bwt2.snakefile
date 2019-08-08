#MODULE: Align fastq files to genome - BOWTIE2 specific calls
#PARAMETERS:
# _logfile=output_path + "/logs/align.log"
#_bwa_q="5"
#_bwa_l="32"
#_bwa_k="2"
_bwt2_threads=8
_samtools_threads=4

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

def bwt2_aln_inputs(inputs):
    """Checks wether the bwt2_aln input is PE or SE
    RETURNS appropriate bwt2 input settings"""
    if len(inputs) > 1:
        #PE
        return ("-1 %s -2 %s" % (inputs[0], inputs[1]))
    else:
        #SE
        return "-U %s" % inputs
                

rule bwt2_aln:
    input:
        getFastq
    output:
        temp(output_path + "/align/{sample}/{sample}.sam")
    params:
        index=config['bwt2_index'],
        _inputs = lambda wildcards, input: bwt2_aln_inputs(input),
        read_group = lambda wildcards: "--rg-id %s --rg \"SM:%s\" --rg \"PL:ILLUMINA\"" % (wildcards.sample, wildcards.sample)
    threads: _bwt2_threads
    message: "ALIGN: Running Bowtie2 alignment"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_bwt2.yaml"
    shell:
        "bowtie2 -p {threads} {params.read_group} -x {params.index} {params._inputs} -S {output} 2>>{log}"

rule bwt2_convert:
    input:
        sam=output_path + "/align/{sample}/{sample}.sam"
    output:
        temp(output_path + "/align/{sample}/{sample}.bam")
    threads: _samtools_threads
    message: "ALIGN: convert sam to bam"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_bwt2.yaml"
    shell:
        "samtools view -@ {threads} -Sb {input} > {output}"

