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
        ls.append("analysis/align/%s/%s_4M_nonChrM.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_4M_unique_nonChrM.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_nonChrM_stat.txt" % (sample,sample))
        ls.append("analysis/frag/%s/%s_fragModel.R" % (sample,sample))
        ls.append("analysis/frag/%s/%s_fragDist.png" % (sample,sample))
    for run in config["runs"].keys():
        ls.append("analysis/frips/%s/%s_frip.txt" % (run,run))
    ls.append("analysis/frips/pbc.csv")
    ls.append("analysis/frips/nonChrM_stats.csv")
    ls.append("analysis/frag/fragSizes.csv")
    ls.append("analysis/frips/frips.csv")
    return ls

def frip_getTreatBam(wildcards):
    """RETURNS the associated 4M_nonChrM.bam for the run's treatment sample"""
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = r[0]
    #print(first)
    ret = "analysis/align/%s/%s_4M_nonChrM.bam" % (first, first)
    return ret

rule frips_all:
    input:
        frips_targets

#SECTION: generate 4M sub-sampled bams
#It's a bit redundant, but I tried my best to reduce the redundancy
rule create_nonChrM:
    """RULE that will generate the base file for both 
    4M_nonChrM.bam AND 4M_unique_nonChrM.bam so we save some redundancy here"""
    input:
        "analysis/align/{sample}/{sample}.sorted.bam"
    params:
        #hack to get the regex in to filter out chrM, random, chrUn
        regex="\'/chrM/d;/random/d;/chrUn/d\'",
    message: "FRiPs: creating the nonChrM SAM file"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_nonChrM.sam'
    shell:
        "samtools view -h {input} | sed -e {params.regex} > {output} 2>>{log}"

rule sample_4M_from_nonChrM:
    """Sample 4M reads from nonChrM SAM file (from create nonChrM)
    ref: https://sourceforge.net/p/samtools/mailman/message/29011091/ 
    see '-s 21.5'
    """
    input:
        'analysis/align/{sample}/{sample}_nonChrM.sam'
    params:
        n="4000000"
    message: "FRiPs: sample- 4M from non-chrM reads"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_4M_nonChrM.bam'
    shell:
        """
        chips/modules/scripts/frips_sample.sh -n {params.n} -i {input} -o {output} 2>>{log}
        """

rule create_unique_nonChrM:
    """Generate _unique_nonChrM.bam by
    Filter out non-uniquely mapped reads from _nonChrM.sam
    """
    input:
        "analysis/align/{sample}/{sample}_nonChrM.sam"
    message: "FRiPs: create uniquely mapped non-chrM reads"
    log:_logfile
    output:
        #make temp
        'analysis/align/{sample}/{sample}_unique_nonChrM.bam'
    shell:
        "samtools view -b -h -F 4 {input} > {output} 2>>{log}"

rule sample_4M_from_uniqueNonChrM:
    """Sample 4M reads from uniqueNonChrM reads
    ref: https://sourceforge.net/p/samtools/mailman/message/29011091/ 
    see '-s 21.5'
    """
    input:
        'analysis/align/{sample}/{sample}_unique_nonChrM.bam'
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

rule frip_calculate:
    """Calculate the frip score"""
    #TODO: if there are more than 1 treatment, merge them??!
    input:
        treat=frip_getTreatBam,
        bed="analysis/peaks/{run}/{run}_sorted_peaks.narrowPeak",
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
        expand("analysis/frips/{sample}/{sample}_pbc.txt", sample=sorted(config["samples"]))
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
        nonChrM_bam="analysis/align/{sample}/{sample}_unique_nonChrM.bam"
    output:
        #temp("analysis/align/{sample}/{sample}_mapping.txt")
        "analysis/align/{sample}/{sample}_nonChrM_stat.txt"
    message: "ALIGN: get nonChrM mapping stats for each bam"
    log: _logfile
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "samtools view -c {input.uniq_bam} > {output} 2>>{log}"
        " && samtools view -c {input.nonChrM_bam} >> {output} 2>> {log}"

rule collect_nonChrM_stats:
    """Aggregate all nonChrM stats for ALL of the samples"""
    input:
        expand("analysis/align/{sample}/{sample}_nonChrM_stat.txt", sample=sorted(config["samples"]))
    output:
        "analysis/frips/nonChrM_stats.csv"
    message: "FRiPs: collect and parse ALL nonChrM stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_collectNonChrM.py -f {files} > {output} 2>>{log}")

rule generate_FragSizeModel:
    """Call macs2 predictd to calculate the fragment size model 
    which will be used to calculate the median fragment sizes for each sample
    see calculate_FragSizes below.

    NOTE: this should be for treatments only but I'm doing this for controls 
    too"""
    input:
        "analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    params:
        #NOTE: this macs2 param can be tricky!
        genome_size = config['genome_size'],
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
    message: "FRiPs: run macs2 fragment size estimation (predictd)"
    log:_logfile
    output:
        #make temp
        "analysis/frag/{samples}/{sample}_fragModel.R"
    shell:
        "{params.pypath} {config[macs2_path]} predictd -i {input} --rfile {output} -g {params.genome_size} 2>>{log}"

rule calculate_FragSizes:
    """Given a macs2 predictd fragment size model (or a set of them)
    calculates the estimated fragment sizes
    """
    input:
        expand("analysis/frag/{sample}/{sample}_fragModel.R", sample=sorted(config['samples']))
    message: "FRiPs: calculate fragment sizes"
    log:_logfile
    output:
        "analysis/frag/fragSizes.csv"
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frag_estFragSize.py -f {files} > {output} 2>>{log}")

rule get_SampleFragLength:
    """Dump all of the sample's fragment lengths into 
    frag/{sample}/{sample}_frags.txt, so we can generate the distribution plot in make_FragPlot
    """
    input:
        "analysis/align/{sample}/{sample}_unique.sorted.bam"
    params:
        awk_cmd = """awk ' $1 <= 1000 && $1 > 0 '"""
    message: "FRAG: get fragment sizes"
    log:_logfile
    output:
        temp("analysis/frag/{sample}/{sample}_frags.txt")
    shell:
        #GRAB out the 9th column, ensuring it's in 1-1000
        "samtools view {input} | cut -f 9 | {params.awk_cmd} > {output} 2>>{log}"
 
rule make_FragPlot:
    """plot the fragment distribution:
    generate the R plot by running frag_plotFragDist.R on _frags.txt
    """
    input:
        "analysis/frag/{sample}/{sample}_frags.txt"
    params:
        name= lambda wildcards: wildcards.sample
    message: "FRAG: plot fragment size distribution plot"
    log:_logfile
    output:
        "analysis/frag/{sample}/{sample}_fragDist.png"
    shell:
        #RUN the R script to get the plot
        "chips/modules/scripts/frag_plotFragDist.R {input} {output} {params.name} 2>>{log}"
    
rule getFripStats:
    """Collect the frips statistics from analysis/frips/{run}/{run}_frip.txt"""
    input:
        expand("analysis/frips/{run}/{run}_frip.txt", run=sorted(config["runs"]))
    output:
        "analysis/frips/frips.csv"
    message: "FRiPs: collecting frips stats for each run"
    log:_logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_getFrips.py -f {files} -o {output} 2>>{log}")
