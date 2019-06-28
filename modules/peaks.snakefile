#MODULE: PEAK CALLING using macs2

#PARAMETERS
# _logfile="analysis/logs/peaks.log"
_macs_fdr="0.01"
_macs_keepdup="1"
_macs_extsize="146"

def getTreats(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])
    #print(rep_n)

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    treat = 2**(rep_n-1) if rep_n > 1 else 0
    r = config['runs'][wildcards.run]
    #print(r)
    #print(treat)
    if treat < len(r) and r[treat]:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[treat],r[treat])]
    #print("TREAT: %s" % tmp)
    if not tmp:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
    return tmp

def getConts(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    cont = 2**(rep_n-1) + 1 if rep_n > 1 else 1
    r = config['runs'][wildcards.run]
    #print(r)
    #print(cont)
    if cont < len(r) and r[cont]:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[cont],r[cont])]
    #print("CONT: %s" % tmp)
    return tmp

def getFilteredTreats(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])
    #print(rep_n)

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    treat = 2**(rep_n-1) if rep_n > 1 else 0
    r = config['runs'][wildcards.run]
    #print(r)
    #print(treat)
    if r[treat]:
        treatSample = config["samples"][r[treat]]
        if len(treatSample) > 1 and ('cutoff' in config) and config['cutoff']:
            #PE
            if treat < len(r):
                tmp = ["analysis/align/%s/%s_unique.sorted.dedup.sub%s.bam" % (r[treat],r[treat],config['cutoff'])]
            #print("TREAT: %s" % tmp)
            if not tmp:
                #NOTE: I can't figure out a proper kill command so I'll do this
                tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
        else:
            #SE
            if treat < len(r):
                tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[treat],r[treat])]
            #print("TREAT: %s" % tmp)
            if not tmp:
                #NOTE: I can't figure out a proper kill command so I'll do this
                tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
    return tmp

def getFilteredConts(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    cont = 2**(rep_n-1) + 1 if rep_n > 1 else 1
    r = config['runs'][wildcards.run]
    #print(r)
    #print(cont)
    if r[cont]:
        contSample = config["samples"][r[cont]]
        if len(contSample) > 1 and ('cutoff' in config) and config['cutoff']:
            #PE
            if cont < len(r):
                tmp = ["analysis/align/%s/%s_unique.sorted.dedup.sub%s.bam" % (r[cont],r[cont],config['cutoff'])]
        else:
            #SE
            if cont < len(r):
                tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[cont],r[cont])]
        #print("CONT: %s" % tmp)
    return tmp

def checkBAMPE(wildcards):
    """Fn returns '-f BAMPE' IF the run's FIRST treatment replicate (sample) is
    Paired-END.
    NOTE: this flag is required for macs2 callpeak, AUTO detect format does not
    work with PE bams!
    """
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = config['samples'][r[0]]
    ret = "-f BAMPE" if len(first) == 2 else ""
    return ret

#NOTE: using the _refs from chips.snakefile
def peaks_targets(wildcards):
    """Generates the targets for this module"""
    #print(wildcards)
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
                ls.append("analysis/peaks/%s/%s_sorted_peaks.broadPeak" % (runRep,runRep))
                ls.append("analysis/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep))
                # ls.append("analysis/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep))
            else:
                ls.append("analysis/peaks/%s/%s_sorted_peaks.narrowPeak" % (runRep,runRep))
                ls.append("analysis/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep))
                ls.append("analysis/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_treat_pileup.bw" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_control_lambda.bw" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_treat_pileup.sorted.bdg.gz" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_control_lambda.sorted.bdg.gz" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_treatment.igv.xml" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_peaks.bed" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_model.R" % (runRep,runRep))
            if ('cutoff' in config) and config['cutoff']:
                if ("macs2_broadpeaks" not in config) or (config["macs2_broadpeaks"] == False):
                    ls.append("analysis/peaks/%s/%s.sub%s_summits.bed" % (runRep,runRep,config["cutoff"]))
                ls.append("analysis/peaks/%s/%s.sub%s_peaks.bed" % (runRep,runRep,config["cutoff"]))
                ls.append("analysis/peaks/%s/%s.sub%s_treat_pileup.bw" % (runRep,runRep,config["cutoff"]))
                ls.append("analysis/peaks/%s/%s.sub%s_control_lambda.bw" % (runRep,runRep,config["cutoff"]))
    ls.append("analysis/peaks/peakStats.csv")
    ls.append("analysis/peaks/run_info.txt")
    ls.append("analysis/peaks/all_treatments.igv.xml")
    return ls

rule peaks_all:
    input:
        peaks_targets

if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
    rule macs2_callpeaks_broad:
        input:
            treat=getTreats,
            cont=getConts
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.broadPeak",
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.gappedPeak",
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.xls",
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}_treat_pileup.bdg"),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}_control_lambda.bdg"),
        params:
            fdr=_macs_fdr,
            keepdup=_macs_keepdup,
            extsize=_macs_extsize,
            genome_size=config['genome_size'],
            outdir="analysis/peaks/{run}.{rep}/",
            name="{run}.{rep}",
            #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
            BAMPE = lambda wildcards: checkBAMPE(wildcards),
            pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
            treatment = lambda wildcards, input: [" -t %s" % i for i in input.treat] if input.treat else "",
            control = lambda wildcards, input: [" -c %s" % i for i in input.cont] if input.cont else "",
        message: "PEAKS: calling peaks with macs2"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
           "{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} "
           "--extsize {params.extsize} --nomodel {params.treatment} {params.control} --broad --broad-cutoff {params.fdr} --outdir {params.outdir} -n {params.name} "

    rule macs2_filtered_callpeaks_broad:
        input:
            treat=getFilteredTreats,
            cont=getFilteredConts
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.broadPeak" % str(config['cutoff']),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.gappedPeak" % str(config['cutoff'])),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.xls" % str(config['cutoff'])),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_treat_pileup.bdg" % str(config['cutoff'])),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_control_lambda.bdg" % str(config['cutoff'])),
        params:
            fdr=_macs_fdr,
            keepdup=_macs_keepdup,
            extsize=_macs_extsize,
            genome_size=config['genome_size'],
            outdir="analysis/peaks/{run}.{rep}/",
            name="{run}.{rep}.sub%s" % str(config['cutoff']),
            #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
            BAMPE = lambda wildcards: checkBAMPE(wildcards),
            pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
            treatment = lambda wildcards, input: [" -t %s" % i for i in input.treat] if input.treat else "",
            control = lambda wildcards, input: [" -c %s" % i for i in input.cont] if input.cont else "",
        message: "PEAKS: calling broad peaks with macs2"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
           "{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} "
           "--extsize {params.extsize} --nomodel {params.treatment} {params.control} --broad --broad-cutoff {params.fdr} --outdir {params.outdir} -n {params.name} "

    rule unsortBoardPeaksToBed:
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.broadPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.bed"
        message: "PEAKS: Converting unsorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule sortBroadPeaks:
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.broadPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.broadPeak"
        message: "PEAKS: sorting the narrowPeaks by -log10qval (col9)"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "sort -r -n -k 9 {input} > {output} "#2>>{log}"

    rule boardPeakToBed:
        """Convert MACS's narrowPeak format, which is BED12 to BED5"""
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.broadPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
        message: "PEAKS: Converting sorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule filteredBoardPeakToBed:
        """Convert MACS's narrowPeak format, which is BED12 to BED5"""
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.broadPeak" % str(config['cutoff'])
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.bed" % str(config['cutoff'])
        message: "PEAKS: Converting sorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule getBroadStats:
        """Counts  number of peaks, # of 10FC, # of 20FC peaks for each sample"""
        input:
            #Generalized INPUT fn defined in chips.snakefile
            _getRepInput("analysis/peaks/$runRep/$runRep_sorted_peaks.broadPeak")
        params:
            files = lambda wildcards, input: [" -f %s" % i for i in input]
        output:
            "analysis/peaks/peakStats.csv"
        message: "PEAKS: collecting peaks stats for each run"
        # log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cidc_chips/modules/scripts/peaks_getPeakStats.py {params.files} -o {output}"# 2>>{log}"

else:
    rule macs2_callpeaks:
        input:
            treat=getTreats,
            cont=getConts
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.narrowPeak",
            "analysis/peaks/{run}.{rep}/{run}.{rep}_summits.bed",
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.xls",
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}_treat_pileup.bdg"),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}_control_lambda.bdg"),
        params:
            fdr=_macs_fdr,
            keepdup=_macs_keepdup,
            # extsize=_macs_extsize,
            genome_size=config['genome_size'],
            outdir="analysis/peaks/{run}.{rep}/",
            name="{run}.{rep}",
            #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
            BAMPE = lambda wildcards: checkBAMPE(wildcards),
            pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
            treatment = lambda wildcards, input: [" -t %s" % i for i in input.treat] if input.treat else "",
            control = lambda wildcards, input: [" -c %s" % i for i in input.cont] if input.cont else "",
        message: "PEAKS: calling peaks with macs2"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
           "{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} "
           "{params.treatment} {params.control} --outdir {params.outdir} -n {params.name} " #2>>{log}"
        #run:
        #    treatment = "-t %s" % input.treat if input.treat else "",
        #    control = "-c %s" % input.cont if input.cont else "",        
        #    shell("{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} --extsize {params.extsize} --nomodel {treatment} {control} --outdir {params.outdir} -n {params.name} 2>>{log}")

    rule macs2_filtered_callpeaks:
        input:
            treat=getFilteredTreats,
            cont=getFilteredConts
        output:
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.narrowPeak" % str(config['cutoff'])),
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_summits.bed" % str(config['cutoff']),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.xls" % str(config['cutoff'])),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_treat_pileup.bdg" % str(config['cutoff'])),
            temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_control_lambda.bdg" % str(config['cutoff'])),
        params:
            fdr=_macs_fdr,
            keepdup=_macs_keepdup,
            # extsize=_macs_extsize,
            genome_size=config['genome_size'],
            outdir="analysis/peaks/{run}.{rep}/",
            name="{run}.{rep}.sub%s" % str(config['cutoff']),
            #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
            BAMPE = lambda wildcards: checkBAMPE(wildcards),
            pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
            treatment = lambda wildcards, input: [" -t %s" % i for i in input.treat] if input.treat else "",
            control = lambda wildcards, input: [" -c %s" % i for i in input.cont] if input.cont else "",
        message: "PEAKS: calling filtered reads peaks with macs2"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
           "{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} "
           "{params.treatment} {params.control} --outdir {params.outdir} -n {params.name} "#2>>{log}"

    rule unsortPeaksToBed:
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.narrowPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.bed"
        message: "PEAKS: Converting unsorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule peakToBed:
        """Convert MACS's narrowPeak format, which is BED12 to BED5"""
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
        message: "PEAKS: Converting sorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule filteredPeakToBed:
        """Convert MACS's narrowPeak format, which is BED12 to BED5"""
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.narrowPeak" % str(config['cutoff'])
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_peaks.bed" % str(config['cutoff'])
        message: "PEAKS: Converting sorted peak file to bed file"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cut -f1,2,3,4,9 {input} > {output} "#2>>{log}"

    rule sortNarrowPeaks:
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.narrowPeak"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak"
        message: "PEAKS: sorting the narrowPeaks by -log10qval (col9)"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "sort -r -n -k 9 {input} > {output} "#2>>{log}"

    rule getPeaksStats:
        """Counts  number of peaks, # of 10FC, # of 20FC peaks for each sample"""
        input:
            #Generalized INPUT fn defined in chips.snakefile
            _getRepInput("analysis/peaks/$runRep/$runRep_sorted_peaks.narrowPeak")
        params:
            files = lambda wildcards, input: [" -f %s" % i for i in input]
        output:
            "analysis/peaks/peakStats.csv"
        message: "PEAKS: collecting peaks stats for each run"
        # log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "cidc_chips/modules/scripts/peaks_getPeakStats.py {params.files} -o {output}"# 2>>{log}"

    rule sortSummits:
        input:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_summits.bed"
        output:
            "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_summits.bed"
        message: "PEAKS: sorting the summits bed by score"
        log:"analysis/logs/peaks/{run}.{rep}.log"
        conda: "../envs/peaks/peaks.yaml"
        shell:
            "sort -r -n -k 5 {input} > {output} "#2>>{log}"



rule macs2_get_fragment:
    input:
        treat=getTreats,
        cont=getConts
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_model.R",
    params:
        #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        treatment = lambda wildcards, input: [" -i %s" % i for i in input.treat] if input.treat else "",
    message: "PEAKS: Get fragment size with macs2"
    log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
       "{params.pypath} {config[macs2_path]} predictd {params.treatment} --rfile {output} -g 'hs' "#2>>{log}"

rule sortBedgraphs:
    """Sort bed graphs--typically useful for converting bdg to bw"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bdg"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg"
    params:
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: sorting bdg pileups {params.msg}"
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "bedSort {input} {output}"# 2>>{log}"

rule bdgToBw:
    """Convert bedGraphs to BigWig"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bw"
    params:
        chroms=config['chrom_lens'],
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: Convert bedGraphs to BigWig {params.msg}"
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "bedGraphToBigWig {input} {params.chroms} {output} "#2>>{log}"

rule sortFilteredBedgraphs:
    """Sort bed graphs--typically useful for converting bdg to bw"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_{suffix}.bdg" % str(config['cutoff'])
    output:
        temp("analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_{suffix}.sorted.bdg" % str(config['cutoff']))
    params:
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: sorting bdg pileups {params.msg}"
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "bedSort {input} {output}"# 2>>{log}"

rule filteredBdgToBw:
    """Convert bedGraphs to BigWig"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_{suffix}.sorted.bdg" % str(config['cutoff'])
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}.sub%s_{suffix}.bw" % str(config['cutoff'])
    params:
        chroms=config['chrom_lens'],
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: Convert bedGraphs to BigWig {params.msg}"
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "bedGraphToBigWig {input} {params.chroms} {output}"# 2>>{log}"

rule gzip_bdg:
    """Space saving rule to compress the bdg output"""
    input:
        bdg="analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg",
        #NOTE: the .bw is NOT used, but it helps ensure rule bdgToBw runs first
        bw="analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bw"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg.gz"
    params:
        #msg just for message below
        msg= lambda wildcards: "%s.%s" % (wildcards.run, wildcards.rep)
    message: "PEAKS: compressing sorted.bdg {params.msg}"
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "gzip {input.bdg} " #2>> {log}"


rule macsRunInfo:
    """Dump the current version of macs and the fdr used into a text file 
    for the report"""
    params:
        fdr = _macs_fdr,
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
    output:
        #MAKE temp
        "analysis/peaks/run_info.txt"
    message: "PEAKS/REPORT - collection macs version and fdr info"
    conda: "../envs/peaks/peaks.yaml"
    shell:
        "{params.pypath} {config[macs2_path]} --version 2> {output} && echo fdr {params.fdr} >> {output}"
    
rule generate_IGV_session:
    """Generates analysis/peaks/all_treatments.igv.xml, a igv session of all
    of the treatment.bw files"""
    input:
        _getRepInput("analysis/peaks/$runRep/$runRep_treat_pileup.bw")
    params:
        genome=config['assembly'],
        treats = lambda wildcards, input: [" -t %s" % i for i in input]
    # log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    output:
        "analysis/peaks/all_treatments.igv.xml"
    message: "PEAKS: generate IGV session for all treatment.bw files"
    shell:
        "cidc_chips/modules/scripts/peaks_generateIGVSession.py -g {params.genome} {params.treats} -o {output}"# 2>>{log}"
        
rule generate_IGV_perTrack:
    """Generates analysis/peaks/{runRep}/{runRep}.igv.xml, a igv session of 
    the treatment.bw file
    **VERY similar to generate_IGV_session, but this is for each individual
    treatment"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_treat_pileup.bw"
    params:
        genome=config['assembly']
    log:"analysis/logs/peaks/{run}.{rep}.log"
    conda: "../envs/peaks/peaks.yaml"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_treatment.igv.xml"
    message: "PEAKS: generate IGV session for {run}.{rep} treatment.bw file"
    shell:
        #NOTE: difference with this call and with generate_IGV_session is we pass the -l param which changes the file path
        "cidc_chips/modules/scripts/peaks_generateIGVSession.py -g {params.genome} -t {input} -o {output} -l"# 2>>{log}"
