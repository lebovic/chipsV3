#MODULE: fastqc- sequence quality scores
#RULES:
#    sample_fastq: subsample 100k reads from each sample to perform fastqc analysis
_logfile="analysis/logs/fastqc.log"

def fastqc_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/fastqc/%s_perSeqGC.txt" % sample)
        ls.append("analysis/fastqc/%s_perSeqQual.txt" % sample)
        ls.append("analysis/fastqc/%s.csv" % sample)
    return ls

def getFastqcInput(wildcards):
    """Get the input file for fastqc. It's either the 100k.fastq sample or the 
     original bam"""
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    if first_file.endswith('.bam'):
        #CLEANER to check for .bam vs (.fastq, fq, fq.gz, etc)
        #ret = [first_file]
        #HACK to get fastqc to name things correctly.  IF we returned
        #the bam file, fastqc will name the output files accordingly
        ret = ["analysis/align/%s/%s_100k.bam" % (s,s)]
    else:
        #HACK: need to give the EVALUATED cannonical path
        ret = ["analysis/align/%s/%s_100k.fastq" % (s,s)]
    #print(ret)
    return ret

def getFastqcBam(wildcards):
    """Used in conjunction w/ getFastqcInput to just return the sample bam
    """
    s = wildcards.sample
    ret = [config["samples"][wildcards.sample][0]]
    #print(ret)
    return ret

def getFastq3(wildcards):
    """Get associated fastqs for each sample.  
    NOTE: if PE, take the first pair"""
    s = config["samples"][wildcards.sample]
    return s[0]

rule fastqc_all:
    input:
        fastqc_targets

rule sample_fastq:
    """Subsample fastq"""
    input:
        getFastq3
    output:
        #MAKE temp
        "analysis/align/{sample}/{sample}_100k.fastq"
    params:
        seed=11,
        #how many to sample
        size=100000
    message: "FASTQC: sample_fastq"
    log:_logfile
    shell:
        "seqtk sample -s {params.seed} {input} {params.size} > {output} 2>>{log}"

rule makeBamLink:
    """USED only when the sample-input is a bam file.  Needed so we name
    the fastqc output files correctly.
    Creates a link to the input bam, and names it {sample}_100k.bam"""
    input:
        getFastqcBam
    output:
        "analysis/align/{sample}/{sample}_100k.bam"
    message: "FASTQC: linking input bam as 100k bam"
    log:_logfile
    shell:
        #NEED to make the intermediary file foo to get it to work
        #DOESN'T work, need abs paths!
        #"ln -s {input} ./foo.bam && mv ./foo.bam {output}"
        "ln -s \"$(readlink -f {input})\" {output}"

    
rule call_fastqc:
    """CALL FASTQC on each sub-sample"""
    input:
        #"analysis/align/{sample}/{sample}_100k.fastq"
        getFastqcInput
    output:
        #"analysis/fastqc/{sample}_100k_fastqc"
        #MAKE temp
        "analysis/fastqc/{sample}_100k_fastqc/fastqc_data.txt",
        #"analysis/fastqc/{sample}_100k_fastqc/
        temp("analysis/fastqc/{sample}_100k_fastqc.html"),
        temp("analysis/fastqc/{sample}_100k_fastqc.zip")
    #threads:
    message: "FASTQC: call fastqc"
    log:_logfile
    shell:
        "fastqc {input} --extract -o analysis/fastqc 2>>{log}"

rule get_PerSequenceQual:
    """extract per sequence quality from fastqc_data.txt"""
    input:
        "analysis/fastqc/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        "analysis/fastqc/{sample}_perSeqQual.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence quality'"
    message: "FASTQC: get_PerSequenceQual"
    log:_logfile
    shell:
        "chips/modules/scripts/fastqc_data_extract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule get_PerSequenceGC:
    """extract per sequence GC contentfrom fastqc_data.txt"""
    input:
        "analysis/fastqc/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        "analysis/fastqc/{sample}_perSeqGC.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence GC content'"
    message: "FASTQC: get_PerSequenceGC"
    log:_logfile
    shell:
        "chips/modules/scripts/fastqc_data_extract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule extract_FastQCStats:
    """extract per sequence GC content, and seq qual stats from fastqc run"""
    input:
        gc = "analysis/fastqc/{sample}_perSeqGC.txt",
        qual = "analysis/fastqc/{sample}_perSeqQual.txt"
    output:
        "analysis/fastqc/{sample}.csv"
    message: "FASTQC: extract_FastQCStats"
    log:_logfile
    shell:
        "chips/modules/scripts/fastqc_stats.py -a {input.qual} -b {input.gc} > {output} 2>>{log}"

