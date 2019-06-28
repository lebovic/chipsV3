import os

def checking_targets(wildcards):
    ls = []
    ls.append("analysis/logs/empty_file_list.txt")
    return ls


def emptyCheckingInput(wildcards):
    # Same as all_targets. Should be modified if all_targets modified.
    # however DO NOT have checking_targets
    _qdnaseq = config["cnv_analysis"]
    ls = []
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(peaks_targets(wildcards))
    ls.extend(fastqc_targets(wildcards))
    ls.extend(conservation_targets(wildcards))
    ls.extend(ceas_targets(wildcards))
    ls.extend(frips_targets(wildcards))
    ls.extend(regulatory_targets(wildcards))
    #Check to see if motif is enabled
    if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] == False:
        if 'motif' in config:
            ls.extend(motif_targets(wildcards))

    #HANDLE CNV/qdnaseq analysis
    if _qdnaseq:
        #ls.extend(qdnaseq_targets(wildcards))

        #check for some inputs
        hasInput = False

        #HACK: for some reason, using the following line causes errors
        #for (run, ls) in config['runs'].items():
        #SO we call getRuns (from above) using a simplified config
        tmp_config = {'metasheet': config['metasheet']}
        runs = getRuns(tmp_config)['runs'].copy()
        for (run) in runs.keys():
            #NOTE: if i do this, this is an error!
            #ls = runs[run]
            if runs[run][1] or runs[run][3]:
                #these are the control sample indices
                hasInput = True
                break
        if hasInput:
            ls.extend(qdnaseq_targets(wildcards))
    # skip running modules that useless in cistrome db 
    if 'ChilinApi' in config and config['ChilinApi'] == True:
        ls.extend(json_targets(wildcards))
        ls.extend(chilin_targets(wildcards))
    else:
        ls.extend(contamination_targets(wildcards))
        ls.extend(mapmaker_targets(wildcards))
        #ls.extend(bam_snapshots_targets(wildcards))
        ls.extend(report_targets(wildcards))
        if "epicypher_analysis" in config and config["epicypher_analysis"]:
            ls.extend(epicypher_targets(wildcards))
    return ls


rule emptyChecking:
    input:
        emptyCheckingInput
    output:
        file="analysis/logs/empty_file_list.txt"
    message:
        "EMPTYCHECKING: checking whether any files are empty"
    run:
        empty_list = []
        for i in input:
            size = os.path.getsize(i)
            if size == 0:
                empty_list.append(i)
            else:
                continue
        with open(output.file,"w") as op:
            op.write("\n".join(empty_list))














