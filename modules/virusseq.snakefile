# Len Taing (2023) TGBTG
#MODULE: virusseq- quantify virusseq viral baits

#_logfile="analysis/logs/virusseq.log"
_threads=8
_viral_bait_chr_name = _vChrName = config.get('virusseq_chrName', "chrVirus")

def virusseq_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #if config.get('virusseq_analysis', False):
    for sample in config["samples"]:
        ls.append(output_path + "/virusseq/%s/%s.virusseq.bam" % (sample, sample))
        ls.append(output_path + "/virusseq/%s/%s.virusseq.%s.bam" % (sample, sample, _viral_bait_chr_name))
        ls.append(output_path + "/virusseq/%s/%s.virusseq.bw" % (sample, sample))
        ls.append(output_path + "/virusseq/%s/%s.virusseq.counts.filtered.txt" % (sample, sample))
    ls.append(output_path + "/virusseq/virusseq.allSample.counts.csv")
    ls.append(output_path + "/virusseq/virusseq.allSample.tpms.csv")
    return ls

rule virusseq_all:
    input:
        virusseq_targets

rule virusseq_getUnmappedReads:
    input:
        output_path + "/align/{sample}/{sample}.sorted.bam"
    output:
        output_path + "/virusseq/{sample}/{sample}.unmapped.bam"
    threads: _threads
    shell:
        "samtools view -b -f 4 -@ {threads} {input} > {output}"

rule virusseq_sortByReadname:
    input:
        output_path + "/virusseq/{sample}/{sample}.unmapped.bam"
    output:
        output_path + "/virusseq/{sample}/{sample}.unmapped.sorted.bam"
    threads: _threads
    shell:
        "sambamba sort -o {output} -N -t {threads} {input}"


rule virusseq_bamToFastq:
    input:
        output_path + "/virusseq/{sample}/{sample}.unmapped.sorted.bam"
    params:
        #fq2 = lambda wildcards: "-fq2 analysis/virusseq/%s/%s.unmapped.fq2" % (wildcards.sample, wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else "",
        fq2 = lambda wildcards: "-2 analysis/virusseq/%s/%s.unmapped.fq2" % (wildcards.sample, wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else "",
    output:
        fq1=output_path + "/virusseq/{sample}/{sample}.unmapped.fq",
        #fq2=output_path + "/virusseq/{sample}/{sample}.unmapped.fq2",
    threads: _threads
    shell:
        #"bamToFastq -i {input} -fq {output.fq1} -fq2 {output.fq2}"
        #"bamToFastq -i {input} -fq {output.fq1} {params.fq2} 2> /dev/null"
        "samtools fastq -@ {threads} -1 {output.fq1} {params.fq2} -0 /dev/null -s /dev/null {input}"

#def virusseq_gzip_input(wildcards):
#    sample = wildcards.sample
#    ls = [output_path + "/virusseq/%s/%s.unmapped.fq" % (sample, sample)]
#    if len(config["samples"][sample]) == 2:
#        ls.append(output_path + "/virusseq/%s/%s.unmapped.fq2" % (sample, sample))
#    return ls

rule virusseq_gzip:
    input:
        #virusseq_gzip_input
        fq1=output_path + "/virusseq/{sample}/{sample}.unmapped.fq",
        #fq2=output_path + "/virusseq/{sample}/{sample}.unmapped.fq2",
    params:
        fq2 = lambda wildcards: output_path + "/virusseq/%s/%s.unmapped.fq2" % (wildcards.sample, wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else "",
    output:
        gz_fq1 = output_path + "/virusseq/{sample}/{sample}.unmapped.fq.gz",
        #DON'T check that the fq2 is also gzipped!
        #gz_fq2 = output_path + "/virusseq/{sample}/{sample}.unmapped.fq2.gz"
    shell:
        #"gzip {input.fq1} {input.fq2}"
        "gzip {input.fq1} {params.fq2}"

#def virusseq_map_input(wildcards):
#    sample = wildcards.sample
#    ls = [output_path + "/virusseq/%s/%s.unmapped.fq.gz" % (sample, sample)]
#    if len(config["samples"][sample]) == 2:
#        ls.append(output_path + "/virusseq/%s/%s.unmapped.fq2.gz" % (sample, sample))
#    return ls

###############################################################################
#ALIGN with BWA
###############################################################################
# rule virusseq_map:
#     input:
#         #virusseq_map_input
#         gz_fq1 = output_path + "/virusseq/{sample}/{sample}.unmapped.fq.gz",
#         #gz_fq2 = output_path + "/virusseq/{sample}/{sample}.unmapped.fq2.gz"
#     params:
#         #bwa_virusseq = "./bwa/hg19Virus/hg19Virus.fa",
#         bwa_virusseq = config['virusseq_bwa_index'],
#         gz_fq2 = lambda wildcards: output_path + "/virusseq/%s/%s.unmapped.fq2.gz" \
# % (wildcards.sample, wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else "",
#         mapq_thresh = config.get('virusseq_mapq','30')
#     threads: _threads
#     output:
#         output_path + "/virusseq/{sample}/{sample}.virusseq.bam"
#     shell:
#         #"bwa mem -t {threads} {params.bwa_virusseq} {input.gz_fq1} {params.gz_fq2} | samtools view -Sb - > {output}"
#         "bwa mem -t {threads} {params.bwa_virusseq} {input.gz_fq1} {params.gz_fq2} | sambamba view -S -f bam -F \"mapping_quality >= {params.mapq_thresh}\" /dev/stdin > {output}"

###############################################################################
#ALIGN with Bowtie2
###############################################################################
rule virusseq_bwt2_map:
    input:
        gz_fq1 = output_path + "/virusseq/{sample}/{sample}.unmapped.fq.gz",
    params:
        bwt2_virusseq = config['virusseq_bowtie2_index'],
        gz_fq2 = lambda wildcards: "-2 analysis/virusseq/%s/%s.unmapped.fq2.gz" % (wildcards.sample, wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else "",
        gz_fq1_flag = lambda wildcards: "-1" if len(config["samples"][wildcards.sample]) == 2 else "-U",
        mapq_thresh = config.get('virusseq_mapq','30'),
        num_mismatch = config.get('virusseq_num_mismatch','1'),
    threads: _threads
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.bam"
    shell:
        #FILTER: mapq >= 30, num mismatches <= 1, properly paired
        "bowtie2 -p {threads} -x {params.bwt2_virusseq} {params.gz_fq1_flag} {input.gz_fq1} {params.gz_fq2} | sambamba view -S -f bam -F \"mapping_quality >= {params.mapq_thresh} and [NM] <= {params.num_mismatch} and proper_pair\" /dev/stdin > {output}"

rule virusseq_dedup:
    input:
        output_path + "/virusseq/{sample}/{sample}.virusseq.bam"
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.dedup.bam"
    threads: _threads
    shell:
        "sambamba markdup -t {threads} --remove-duplicates {input} {output}"

rule virusseq_sortByName:
    input:
        #output_path + "/virusseq/{sample}/{sample}.virusseq.bam"
        output_path + "/virusseq/{sample}/{sample}.virusseq.dedup.bam"
        #TODO: if using dedup, then sorted.bam becomes sorted.dedup.bam
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam"
    threads: _threads
    shell:
        "sambamba sort -o {output} -t {threads} {input}"

rule virusseq_indexSorted:
    input:
        output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam"
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam.bai"
    shell:
        "sambamba index {input}"

rule virusseq_extractChrVirus:
    input:
        bam = output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam",
        bai = output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam.bai"
    params:
        chrName = _viral_bait_chr_name,
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.%s.bam" % _viral_bait_chr_name
    shell:
        "samtools view -hb {input.bam} {params.chrName} > {output}"

rule virusseq_sortAndIndexChrVirus:
    input:
        output_path + "/virusseq/{sample}/{sample}.virusseq.%s.bam" % _viral_bait_chr_name
    output:
        bam = output_path + "/virusseq/{sample}/{sample}.virusseq.%s.sorted.bam" % _viral_bait_chr_name,
        bai = output_path + "/virusseq/{sample}/{sample}.virusseq.%s.sorted.bam.bai" % _viral_bait_chr_name,
    threads: _threads        
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input} && samtools index {output.bam}"

#KEY rule: given a .bed file of viral bait regions, this rule will give a
#read count for each virus region/"gene"
rule virusseq_countViralReads:
    input:
        bam = output_path + "/virusseq/{sample}/{sample}.virusseq.%s.sorted.bam" % _viral_bait_chr_name,
        bai = output_path + "/virusseq/{sample}/{sample}.virusseq.%s.sorted.bam.bai" % _viral_bait_chr_name,
    params:
        #four column bed file with chrom, start, end, and virus name for each
        #virus in the bait set
        #viral_bed = "/data/static_libraries/viper/ref_files/hg19/virusseq/hg19Virus.bed",
        viral_bed = config['virusseq_bed'],
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.counts.txt"
    shell:
        "bedtools multicov -bams {input.bam} -bed {params.viral_bed} > {output}"

# LEN: obsolete
# rule virusseq_filteredCounts:
#     """Remove all of the entries in virusseq.counts.txt that have 0 counts"""
#     input:
#         output_path + "/virusseq/{sample}/{sample}.virusseq.counts.txt"
#     params:
#         awk_cmd = "awk \'$5 != 0\'"
#     output:
#         output_path + "/virusseq/{sample}/{sample}.virusseq.counts.filtered.txt"
#     shell:
#         #NOTE: cutting the last two cols from the 5 col bed b/c those are most relevant
#         #NO cutting
#         #"{params.awk_cmd} {input} > {output}"
#         #cutting and sorting by read-count (descending
#         "{params.awk_cmd} {input} | cut -f 4,5 | sort -k 2nr > {output}"

#HACK we're just going to avoid getting uniquely mapped reads by using what we #have
# rule virusseq_mapping_deleteMe:
#     input:
#         output_path + "/align/{sample}/{sample}.sorted.bam"
#     output:
#         output_path + "/virusseq/{sample}/{sample}_mapping.txt",
#     threads: _threads
#     shell:
#         "samtools flagstat -@ {threads} {input} > {output}"
    
rule virusseq_filteredCounts:
    """Remove all of the entries in virusseq.counts.txt that have 0 counts
    and normalize values using TPM
    """
    input:
        counts= output_path +"/virusseq/{sample}/{sample}.virusseq.counts.txt",
        mapping = output_path + "/align/{sample}/{sample}_mapping.txt",
        #HACK input below
        #mapping = output_path + "/virusseq/{sample}/{sample}_mapping.txt",
    output:
        output_path +"/virusseq/{sample}/{sample}.virusseq.counts.filtered.txt"
    shell:
        "cidc_chips/modules/scripts/virusseq_filterAndNormalize.py -f {input.counts} -m {input.mapping} -o {output}"

rule virusseq_aggregateFilteredCts:
    input:
        expand(output_path + "/virusseq/{sample}/{sample}.virusseq.counts.filtered.txt", sample = config['samples'])
    params:
        files = lambda wildcards,input: " -f ".join(input),
        samples =lambda wildcards: " -s ".join([s for s in config['samples']]),
    output:
        cts = output_path + "/virusseq/virusseq.allSample.counts.csv",
        tpm = output_path + "/virusseq/virusseq.allSample.tpms.csv",
    shell:
        "cidc_chips/modules/scripts/virusseq_aggregateFiltered.py -f {params.files} -s {params.samples} -o {output.cts} -t {output.tpm}"

rule virusseq_bamTobedgraph:
    """Create a bedgraph for the chrM reads"""
    input:
        #LEN: moving this from virusseq.chrM.sorted.bam to
        #virusseq.sorted.bam IN case samples don't have any viral contam.
        #output_path + "/virusseq/{sample}/{sample}.virusseq.chrM.sorted.bam"
        output_path + "/virusseq/{sample}/{sample}.virusseq.sorted.bam"
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.bdg"
    shell:
        #LEN: does not sort lexographically, which is what bdgTobw needs
        #"bedtools genomecov -ibam {input} -bg > {output}"
        "bedtools genomecov -ibam {input} -bg | sort -k1,1 -k2,2n > {output}"

rule virusseq_bdgTobw:
    input:
        output_path + "/virusseq/{sample}/{sample}.virusseq.bdg"
    params:
        virusseq_len=config['virusseq_chrom_len'],
    output:
        output_path + "/virusseq/{sample}/{sample}.virusseq.bw"
    shell:
        #"bedGraphToBigWig {input} {params.virusseq_len} {output}"
        "cidc_chips/modules/scripts/virusseq_bdgToBw.sh {input} {params.virusseq_len} {output}"
        
