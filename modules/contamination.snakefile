#MODULE: contamination module to check for sample contamination
_logfile="analysis/logs/contamination.log"
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"
_bwa_threads=4

###############################################################################
# HELPERS
###############################################################################
_contaminationPanel= config['contamination_panel'] if 'contamination_panel' in config and config['contamination_panel'] else []

def extractIndexName(path):
    """Given a contamination panel path, e.g. /some/path/to/BWA/hg19.fa,
    this fn returns the index filename, e.g. hg19"""
    tmp = path.split("/")[-1]
    #lop off the suffix
    return ".".join(tmp.split(".")[:-1])
_contaminationNames = list(map(extractIndexName, _contaminationPanel))
_contaminationDict = dict(zip(_contaminationNames, _contaminationPanel))
###############################################################################

rule contamination_all:
    input:
        expand("analysis/contam/{sample}/{sample}.{panel}.sai", sample=config['samples'].keys(), panel=_contaminationNames),
        expand("analysis/contam/{sample}/{sample}.{panel}.bam", sample=config['samples'].keys(), panel=_contaminationNames),
        expand("analysis/contam/{sample}/{sample}.{panel}.txt", sample=config['samples'].keys(), panel=_contaminationNames),
        expand("analysis/contam/{sample}/{sample}_contamination.txt", sample=config['samples'].keys()),


rule contamination:
    """For each sample, run an alignment for each entry in the 
    contaminationPanel and get the mapping rate"""
    input:
        #from FASTQC: sample_fastq
        fastq="analysis/align/{sample}/{sample}_100k.fastq",
    output:
        #TEMP file
        "analysis/contam/{sample}/{sample}.{panel}.sai"
    params:
        index=lambda wildcards: _contaminationDict[wildcards.panel],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel
    threads: _bwa_threads
    message: "CONTAMINATION: checking {params.sample} against {params.panel}"
    log: _logfile
    shell:
        "bwa aln -q {params.bwa_q} -l {params.bwa_l} -k {params.bwa_k} -t {threads} {params.index} {input.fastq} > {output} 2>>{log}"

rule contaminationToBam:
    """Convert the sai to bams"""
    input:
        sai="analysis/contam/{sample}/{sample}.{panel}.sai",
        fastq="analysis/align/{sample}/{sample}_100k.fastq"
    output:
        #TEMP file
        "analysis/contam/{sample}/{sample}.{panel}.bam"
    params:
        index=lambda wildcards: _contaminationDict[wildcards.panel],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -"
    threads: _bwa_threads
    message: "CONTAMINATION: Converting BWA alignments (.sai) to BAM"
    log: _logfile
    shell:
        #NOTE: all contamination alignments are SE!
        "bwa samse {params.index} {input.sai} {input.fastq} 2>>{log} | samtools {params.hack} > {output} 2>>{log}"

rule contaminationStats:
    """Extract the mapping stats for each species"""
    input:
        "analysis/contam/{sample}/{sample}.{panel}.bam"
    output:
        #TEMP
        "analysis/contam/{sample}/{sample}.{panel}.txt"
    params:
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel,
        #READ out the 5th row, the first element (and divide by 100/100000)
        awk_cmd = "\'NR==5 {print $1 / 1000}\'"
    message: "CONTAMINATION: get mapping stats {params.sample}:{params.panel}"
    log: _logfile
    shell:
        "samtools flagstat {input} | awk {params.awk_cmd} > {output} 2>>{log}"

rule contaminationCollectStats:
    """Collect the mapping stats across the entire panel"""
    input:
        expand("analysis/contam/{{sample}}/{{sample}}.{panel}.txt", panel=_contaminationNames)
    output:
        "analysis/contam/{sample}/{sample}_contamination.txt"
    run:
        for (n, f) in zip(_contaminationNames, input):
            shell("per=$(cat {f}) && echo {n} $per >> {output}")
