#2023-10-17 MODULE originally written by Ruiyang He
#MODULE: Align fastq files to genome - Chromap specific calls
#PARAMETERS:
#_logfile="analysis/logs/align.log"
_chromap_threads=8

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule chromap:
    input:
        getFastq
    output:
        bam=temp(output_path + "/align/{sample}/{sample}.bam"),
        sam=temp(output_path + "/align/{sample}/{sample}.sam"),
    params:
        index=config['chromap_index'],
        ref=config['chromap_ref']
    threads: _chromap_threads
    message: "ALIGN: Running Chromap for alignment"
    log: output_path + "/logs/align_chromap_{sample}.log"
    run:
        # single-end
        if len(input) > 1:
            _inputs=("-1 %s -2 %s" % (input[0], input[1]))
        else:
            _inputs="-1 %s" % input
        shell("chromap -x {params.index} -r {params.ref} {_inputs} -t {_chromap_threads} --SAM -o {output.sam} --trim-adapters -l 2000 -q 0 2>>{log}")
        shell("samtools view -Sb {output.sam} > {output.bam} 2>>{log}")

rule align_chromapInfo:
    """Dump the current version of bwa into a text file for the report"""
    output:
        output_path + "/align/run_info.txt"
    message: "ALIGN/REPORT - collect chromap version info"
    #conda: "../envs/align/align_bwa.yaml"
    shell:
        src_path + "/modules/scripts/align_parseChromapVersion.py -o {output}"
