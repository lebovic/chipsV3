#MODULE: Align fastq files to genome - common rules
_align_threads=4

rule align_all:
    input:
        expand("analysis/align/{sample}/{sample}.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}.sorted.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}_unique.sorted.bam", sample=config["samples"]),
        expand("analysis/align/{sample}/{sample}.unmapped.fq.gz", sample=config["samples"]),
        "analysis/align/mapping.csv",

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
    threads: _align_threads
    shell:
        "samtools sort {input} -o {output} --threads {threads} 2>>{log}"

rule extractUnmapped:
    """Extract the unmapped reads and save as {sample}.unmapped.bam"""
    input:
        "analysis/align/{sample}/{sample}.bam"
    output:
        #THIS should be TEMP
        "analysis/align/{sample}/{sample}.unmapped.bam"
    message: "ALIGN: extract unmapped reads"
    log: _logfile
    threads: _align_threads
    shell:
        #THIS extracts all unmapped reads
        #"samtools view -b -f 4 --threads {threads} {input} >{output} 2>>{log}"
        #THIS extracts all READ (pairs) where at least one in unmapped
        #ref: https://www.biostars.org/p/56246/ search "rgiannico"
        "samtools view -b -F 2 --threads {threads} {input} > {output} 2>>{log}"

rule bamToFastq:
    """Convert the unmapped.bam to fastq"""
    input:
        "analysis/align/{sample}/{sample}.unmapped.bam"
    output:
        "analysis/align/{sample}/{sample}.unmapped.fq"
    params:
        #handle PE alignments!
        mate2 = lambda wildcards: "-fq2 analysis/align/{sample}/{sample}.unmapped.fq2" if len(config["samples"][wildcards.sample]) == 2 else ""
    message: "ALIGN: convert unmapped bam to fastq"
    log: _logfile
    shell:
        "bamToFastq -i {input} -fq {output} {params.mate2}"

rule gzipUnmappedFq:
    """gzip unmapped fq(s)"""
    input:
        "analysis/align/{sample}/{sample}.unmapped.fq"
    output:
        "analysis/align/{sample}/{sample}.unmapped.fq.gz"
    params:
        #handle PE alignments!
        mate2 = lambda wildcards: "analysis/align/{sample}/{sample}.unmapped.fq2" if len(config["samples"][wildcards.sample]) == 2 else ""
    message: "ALIGN: gzip unmapped fq files"
    log: _logfile
    shell:
        "gzip {input} {params} 2>>{log}"
