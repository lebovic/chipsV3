#MODULE: FRiPs- Calculating the fraction of reads under peak (score/percentage)

#helper fns to get samples associated for each run
#copied directly from peaks.snakefile
#TODO: centralize these helper fns!

#PARAMETERS:
_logfile="analysis/logs/frips.log"
_macs_fdr="0.01"
_macs_keepdup="1"
_macs_extsize="146"
_macs_species="hs"

def getTreats(wildcards):
    r = config['runs'][wildcards.run]
    #print(r)
    #convert SAMPLE names to BAMS
    tmp=["analysis/align/%s/%s_4M_unique_nonChrM.bam" % (s,s) for s in r[:2] if s]
    #print(tmp)
    return tmp

def getConts(wildcards):
    r = config['runs'][wildcards.run]
    #print(r)
    #convert SAMPLE names to BAMS
    tmp=["analysis/align/%s/%s_4M_unique_nonChrM.bam" % (s,s) for s in r[2:4] if s]
    #print(tmp)
    return tmp


rule frips_all:
    input:
        expand("analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam", sample=config["samples"]),
        expand("analysis/peaks/{run}/{run}_4M_peaks.narrowPeak", run=config["runs"].keys()),
        expand("analysis/frips/{run}/{run}_frip.txt",run=config["runs"].keys()),

rule sample_unique_nonChrM:
    """Sample uniquely mapped, nonChrM reads from the SAMPLE
    ref: https://www.biostars.org/p/56246/ #for uniquely mapped reads
    ref: https://www.biostars.org/p/128967/ #for the chrM, chrRandom filter
    """
    input:
        "analysis/align/{sample}/{sample}.bam"
    params:
        #hack to get the regex in to filter out chrM, random, chrUn
        regex="\'/chrM/d;/random/d;/chrUn/d\'",
    message: "FRiPs: sample- uniquely mapped non-chrM reads"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_unique_nonChrM.sam'
    shell:
        "samtools view -b -F 4 {input} | sed -e {params.regex} > {output} 2>>{log}"

rule sample_4M_from_uniqueNonChrM:
    """Sample 4M reads from uniqueNonChrM reads
    ref: https://sourceforge.net/p/samtools/mailman/message/29011091/ 
    see '-s 21.5'
    """
    input:
        'analysis/align/{sample}/{sample}_unique_nonChrM.sam'
    params:
        n="4000000"
    message: "FRiPs: sample- 4M from uniquely mapped non-chrM reads"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam'
    shell:
        """
        chips/modules/scripts/frips_sample.sh -n {params.n} -i {input} -o {output} 2>>{log}
        """

rule frips_callpeaks:
    """CALL PEAKS on the 4M reads so we can calc frips"""
    #just handle narrowPeaks for now
    input:
        treat=getTreats,
        cont=getConts
    output:
        "analysis/peaks/{run}/{run}_4M_peaks.narrowPeak",
        "analysis/peaks/{run}/{run}_4M_peaks.xls",
        "analysis/peaks/{run}/{run}_4M_summits.bed",
    params:
        fdr=_macs_fdr,
        keepdup=_macs_keepdup,
        extsize=_macs_extsize,
        species=_macs_species,
        outdir="analysis/peaks/{run}/",
        name="{run}_4M"
    message: "FRiPs: call peaks from sub-sample"
    log:_logfile
    run:
        #NOTE: here's the broadPeak call
        #macs2 callpeak -q [fdr=0.01] --keep-dup [keep_dup=1] --broad -g {param[species]} -t [4M.bam] -c [4M.bam optional] -n [description="SAMPLE_4M"]
        #DIFF: in narrow, --shiftsize and --nomodel is replace by --broad!
        #AND output named broadPeak vs. narrowPeak
        #JOIN treatment and control replicate samples
        treatment = "-t %s" % " ".join(input.treat) if input.treat else "",
        control = "-c %s" % " ".join(input.cont) if input.cont else ""
        shell("macs2 callpeak -q {params.fdr} --keep-dup {params.keepdup} --extsize {params.extsize} --nomodel -g {params.species} {treatment} {control} --outdir {params.outdir} -n {params.name} 2>>{log}")

rule frip_calculate:
    """Calculate the frip score"""
    #TODO: if there are more than 1 treatment, merge them??!
    input:
        treat=getTreats,
        bed="analysis/peaks/{run}/{run}_4M_peaks.narrowPeak",
    output:
        "analysis/frips/{run}/{run}_frip.txt"
    params:
        pval="1E-9"
    message: "FRiPs: calculate frips"
    log:_logfile
    shell:
        "chips/modules/scripts/frips_calculate.sh -a {input.treat} -b {input.bed} -p {params.pval} > {output} 2>>{log}"
