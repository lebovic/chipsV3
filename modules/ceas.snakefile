#MODULE: CEAS- annotating where the peaks fall (promoter, exon, intron, interg)
_logfile="analysis/logs/ceas.log"

def ceas_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        ls.append("analysis/ceas/%s/%s_summary.txt" % (run,run))
        ls.append("analysis/ceas/%s/%s_DHS_peaks.bed" % (run,run))
        ls.append("analysis/ceas/%s/%s_DHS_stats.txt" % (run,run))
        ls.append("analysis/ceas/%s/%s_velcro_peaks.bed" % (run,run))
        ls.append("analysis/ceas/%s/%s_velcro_stats.txt" % (run,run))

    #ADD bam_regionStats
    for sample in config["samples"]:
        if config['exons']:
            ls.append("analysis/ceas/samples/%s/%s.exons" % (sample,sample))
        if config['promoters']:
            ls.append("analysis/ceas/samples/%s/%s.promoters" % (sample,sample))
        if config['DHS']:
            ls.append("analysis/ceas/samples/%s/%s.DHS" % (sample,sample))
    ls.append("analysis/ceas/samples/bamRegionStats.csv")
    ls.append("analysis/ceas/dhs.csv")
    return ls

rule ceas_all:
    input:
        ceas_targets

rule ceas:
    """Annotate peak regions"""
    input:
        "analysis/peaks/{run}/{run}_summits.bed"
    output:
        promoter="analysis/ceas/{run}/{run}_summits_promoter.bed",
        exon="analysis/ceas/{run}/{run}_summits_exon.bed",
        intron="analysis/ceas/{run}/{run}_summits_intron.bed",
        intergenic="analysis/ceas/{run}/{run}_summits_intergenic.bed",
        summary="analysis/ceas/{run}/{run}_summary.txt",
    message: "CEAS: annotating peak regions"
    log: _logfile
    params:
        db=config['geneTable'],
        path="analysis/ceas/{run}/"
    shell:
        "chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} > {output.summary} 2>>{log}"


rule takeTop5k:
    """Take the top 5000 sites"""
    input:
        "analysis/peaks/{run}/{run}_sorted_peaks.narrowPeak.bed"
    params:
        n=5000
    message: "DHS: Take top sites"
    log: _logfile
    output:
        temp('analysis/ceas/{run}/{run}_sorted_5k_peaks.bed')
    shell:
        "head -n {params.n} {input} > {output} 2>>{log}"

#------------------------------------------------------------------------------
rule DHS_intersectDHS:
    """Intersect PEAKS with DHS regions"""
    input:
        'analysis/ceas/{run}/{run}_sorted_5k_peaks.bed'
    params:
        dhs=config['DHS']
    message: "DHS: intersect PEAKS with DHS regions"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_DHS_peaks.bed'
    shell:
        "intersectBed -wa -u -a {input} -b {params.dhs} > {output} 2>>{log}"

rule DHS_stat:
    """collect DHS stats"""
    input:
        n='analysis/ceas/{run}/{run}_sorted_5k_peaks.bed',
        dhs='analysis/ceas/{run}/{run}_DHS_peaks.bed'
    message: "DHS: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_DHS_stats.txt'
    shell:
        "wc -l {input.n} {input.dhs} > {output} 2>>{log}"
#------------------------------------------------------------------------------
rule VELCRO_intersectVelcro:
    """Intersect PEAKS with velcro regions"""
    input:
        'analysis/ceas/{run}/{run}_sorted_5k_peaks.bed'
    params:
        velcro=config['velcro_regions']
    message: "VELCRO: intersect PEAKS with velcro regions"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_velcro_peaks.bed'
    shell:
        "intersectBed -wa -u -a {input} -b {params.velcro} > {output} 2>>{log}"
    
rule VELCRO_stat:
    """collect VELCRO stats"""
    input:
        n='analysis/ceas/{run}/{run}_sorted_5k_peaks.bed',
        velcro='analysis/ceas/{run}/{run}_velcro_peaks.bed'
    message: "VELCRO: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}/{run}_velcro_stats.txt'
    shell:
        "wc -l {input.n} {input.velcro} > {output} 2>>{log}"

rule bam_regionStat:
    """count the number of reads in promoter, exon, dhs--these regions
    are defined in the config.yaml"""
    input:
        "analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    params:
        bed = lambda wildcards: config[wildcards.region]
    message: "CEAS: bam stat region"
    log: _logfile
    output:
        #make temp
        'analysis/ceas/samples/{sample}/{sample}.{region}'
    shell:
        "chips/modules/scripts/meta_bamRegionCount.sh -i {input} -b {params.bed} -o {output} 2>> {log}"

rule collect_BamRegionStats:
    """collect the BAM region stats into a single file"""
    input:
        #INPUT the stats directories--
        #hack just add all of the file so we dont get a missing input exception
        #and collect the directories down below
        dhs = expand("analysis/ceas/samples/{sample}/{sample}.DHS", sample=config['samples']),
        prom = expand("analysis/ceas/samples/{sample}/{sample}.promoters", sample=config['samples']),
        exon = expand("analysis/ceas/samples/{sample}/{sample}.exons", sample=config['samples'])
    message: "CEAS: collect bam region stats"
    log: _logfile
    output:
        'analysis/ceas/samples/bamRegionStats.csv'
    run:
        #REMOVE the file to get directories - 
        #Aribitrarily USE only the exon list
        dirs = ["/".join(d.split("/")[:-1]) for d in input.exon]
        d = " -d ".join(dirs)
        shell("chips/modules/scripts/ceas_collectBamRegStats.py -d {d} > {output} 2>>{log}")

rule collect_DHSstats:
    """collect the DHS stats into a single file"""
    input:
        expand("analysis/ceas/{run}/{run}_DHS_stats.txt", run=config["runs"])
    message: "CEAS: collect DHS stats"
    log: _logfile
    output:
        'analysis/ceas/dhs.csv'
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/peaks_getDHSstats.py -f {files} -o {output} 2>>{log}")
