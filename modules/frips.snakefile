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

def frips_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s_4M_unique_nonChrM.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_nonChrM_stat.txt" % (sample,sample))
    for run in config["runs"].keys():
        ls.append("analysis/peaks/%s/%s_4M_peaks.narrowPeak" % (run,run))
        ls.append("analysis/frips/%s/%s_frip.txt" % (run,run))
    ls.append("analysis/frips/pbc.csv")
    ls.append("analysis/frips/nonChrM_stats.csv")
    return ls

def getTreats(wildcards):
    r = config['runs'][wildcards.run]
    #print(r[:2])
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
        frips_targets

rule sample_unique_nonChrM:
    """Sample uniquely mapped, nonChrM reads from the SAMPLE
    ref: https://www.biostars.org/p/56246/ #for uniquely mapped reads
    ref: https://www.biostars.org/p/128967/ #for the chrM, chrRandom filter
    """
    input:
        #"analysis/align/{sample}/{sample}.bam"
        "analysis/align/{sample}/{sample}_unique.bam"
    params:
        #hack to get the regex in to filter out chrM, random, chrUn
        regex="\'/chrM/d;/random/d;/chrUn/d\'",
    message: "FRiPs: sample- uniquely mapped non-chrM reads"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_unique_nonChrM.sam'
    shell:
        "samtools view -h -F 4 {input} | sed -e {params.regex} > {output} 2>>{log}"

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

rule frip_pbc:
    """Generate the PBC histogram for each normalized sample, which will be 
    used to calculate N1, Nd, and PBC (for the report)
    """
    input:
        "analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    output:
        #make temp
        "analysis/frips/{sample}/{sample}_pbc.txt"
    message: "FRiP: generate PBC histogram for each sample/bam"
    log: _logfile
    shell:
        "chips/modules/scripts/frips_pbc.sh -i {input} -o {output} 2>> {log}"

rule collect_pbc:
    """Collect and parse out the PBC for the ALL of the samples"""
    input:
        expand("analysis/frips/{sample}/{sample}_pbc.txt", sample=config["samples"])
    output:
        "analysis/frips/pbc.csv"
    message: "ALIGN: collect and parse ALL pbc stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_collectPBC.py -f {files} > {output} 2>>{log}")

rule nonChrM_stats:
    """Get the nonChrM mapping stats for each aligment run"""
    input:
        #NOTE: uniq_bam is generated in align_common module-
        #HACK- taking advantage that this moulde is loaded AFTER align_common
        uniq_bam="analysis/align/{sample}/{sample}_unique.bam",
        nonChrM_sam="analysis/align/{sample}/{sample}_unique_nonChrM.sam"
    output:
        #temp("analysis/align/{sample}/{sample}_mapping.txt")
        "analysis/align/{sample}/{sample}_nonChrM_stat.txt"
    message: "ALIGN: get nonChrM mapping stats for each bam"
    log: _logfile
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "samtools view -c {input.uniq_bam} > {output} 2>>{log}"
        " && samtools view -c {input.nonChrM_sam} >> {output} 2>> {log}"

rule collect_nonChrM_stats:
    """Aggregate all nonChrM stats for ALL of the samples"""
    input:
        expand("analysis/align/{sample}/{sample}_nonChrM_stat.txt", sample=config["samples"])
    output:
        "analysis/frips/nonChrM_stats.csv"
    message: "FRiPs: collect and parse ALL nonChrM stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_collectNonChrM.py -f {files} > {output} 2>>{log}")
