#MODULE: json- generating the analysis/json folder
_logfile="analysis/logs/json.log"

def json_targets(wildcards):
    ls = []
    # ls.append("analysis/json")
    for sample in config["samples"]:
        #conserv is run
        ls.append("analysis/json/%s/%s_conserv.json" % (sample, sample))
        
        #contam is sample-- see chilin.snakefile- fastqc as example
        ls.append("analysis/json/%s/%s_contam.json" % (sample, sample))
        
        #DHS is run
        if config["DHS"]:
            ls.append("analysis/json/%s/%s_dhs.json" % (sample, sample))
            
        #Velcro is run
        # if config["velcro_regions"]:
        #     ls.append("analysis/json/%s/%s_velcro.json" % (sample, sample))

        #Not sure about enrich_meta
        ls.append("analysis/json/%s/%s_enrich_meta.json" % (sample, sample))

        #Fastqc is sample-- see chilin.snakefile- fastqc as example
        ls.append("analysis/json/%s/%s_fastqc.json" % (sample, sample))
        #frag is sample-- see chilin.snakefile- fastqc as example
        ls.append("analysis/json/%s/%s_frag.json" % (sample, sample))
        #frip is run
        ls.append("analysis/json/%s/%s_frip.json" % (sample, sample))
        #map is sample---- see chilin.snakefile- fastqc as example
        ls.append("analysis/json/%s/%s_map.json" % (sample, sample))
        # ls.append("analysis/json/%s/%s_macs2_rep.json" % (sample, sample))
        #macs is run
        ls.append("analysis/json/%s/%s_macs2.json" % (sample, sample))
        #meta is run
        ls.append("analysis/json/%s/%s_meta.json" % (sample, sample))
        #pbc is run
        ls.append("analysis/json/%s/%s_pbc.json"% (sample, sample))
        # ls.append("analysis/json/%s/%s_rep.json" % (sample, sample))
    return ls

def json_getRunAndRep(wildcards):
    item = 0
    replicate = []
    for i in range(2,len(r)+2,2):
        a = int(i/2)
        replicate.append("rep%s" % str(a))
        replicate.append("rep%s" % str(a))
        item += 2
    # print(replicate)
    sample = wildcards.sample
    # print(sample)
    key_list=[]
    value_list=[]
    for key,value in config["runs"].items():
        key_list.append(key)
        value_list.append(value)
    for i in value_list:
        if sample in i:
            run = key_list[value_list.index(i)]
            rep = replicate[i.index(sample)]
            break
        else:
            continue
    return "%s.%s" % (run, rep)

# def velcro_target(wildcards):
#     input = "analysis/ceas/%s/%s_velcro_summary.txt" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
#     return input

def enrich_meta_target(wildcards):
    sample = wildcards.sample
    input = []
    meta="analysis/ceas/%s/%s_summary.txt" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    mapped="analysis/align/%s/%s_mapping.txt" % (sample,sample)
    dhs="analysis/ceas/%s/%s_DHS_summary.dhs" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    input.append(meta)
    input.append(mapped)
    input.append(dhs)
    # print(input)
    return input

rule json_all:
    input:
        json_targets

rule json_conservation:
    input:
        lambda wildcards:"analysis/conserv/%s/%s_conserv.txt" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
        #"analysis/conserv/{run}.{rep}/{run}.{rep}_conserv.txt"
    output:
        "analysis/json/{sample}/{sample}_conserv.json"
    params:
        basics = "-b %s " % config["basics"] if "basics" in config else "",
        factor = "-f %s " % config["factor"] if "factor" in config else "",
        TF = "-T %s " % config["TF"] if "TF" in config else ""
    message: "JSON: generate conservation json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_conserv.py -i {input} -o {output} {params.basics}{params.factor}{params.TF}-I {wildcards.sample} 2>>{log}"

rule json_comtamination:
    input:
        "analysis/contam/{sample}/{sample}_contamination.txt"
    output:
        "analysis/json/{sample}/{sample}_contam.json"
    params:

    message: "JSON: generate comtamination json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_comtamination.py -i {input} -o {output} -I {wildcards.sample} -s {wildcards.sample}"


rule json_dhs:
    input:
        lambda wildcards:"analysis/ceas/%s/%s_DHS_summary.dhs" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    output:
        "analysis/json/{sample}/{sample}_dhs.json"
    message: "JSON: generate DHS json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_dhs.py -i {input} -o {output} -I {wildcards.sample}"

# rule json_velcro:
#     input:
#         velcro_target
#     output:
#         "analysis/json/{sample}/{sample}_velcro.json"
#     params:

#     message: "JSON: generate velcro json"
#     log: _logfile
#     shell:
#         "cidc_chips/modules/scripts/json/json_velcro.py  "

rule json_enrich_meta:
    input:
        enrich_meta_target
    output:
        "analysis/json/{sample}/{sample}_enrich_meta.json"
    params:
        has_dhs = "-H %s " % config["DHS"],
        down = True
    message: "JSON: generate meta enrichment json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_enrich_meta.py -i {input[0]} -m {input[1]} -o {output} {params.has_dhs} -D {input[2]} -d {params.down} -I {wildcards.sample} -s {wildcards.sample}"

rule json_fastqc:
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        "analysis/json/{sample}/{sample}_fastqc.json"
    params:
        ids = "-s %s" % config['sample'] if "sample" in config else "",
    message: "JSON: generate fastqc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_fastqc.py -i {input} -o {output} {params.ids} -I {run}"

rule json_frag:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_model.R" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    output:
        json = "analysis/json/{sample}/{sample}_frag.json",
    params:
        sd_R = lambda wildcards: "analysis/peaks/%s/%s_model_sd.R" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    message: "JSON: generate frag json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frag.py -r {input} -o {output} -R {params.sd_R} -f BAMSE -s {wildcards.sample}"

rule json_frip:
    input:
        lambda wildcards: "analysis/frips/%s/%s_frip.txt" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    output:
        "analysis/json/{sample}/{sample}_frip.json"
    message: "JSON: generate frip json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frip.py -i {input} -o {output} -s {wildcards.sample}"

rule json_macs2:
    input:
        lambda wildcards: "analysis/peaks/%s/%s_peaks.xls" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    output:
        "analysis/json/{sample}/{sample}_macs2.json"
    message: "JSON: generate macs2 json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_macs2.py -i {input} -o {output} -I {wildcards.sample}"

# rule json_macs2_rep:
#     input:
#         ""
#     output:
#         "analysis/json/{sample}/{sample}_macs2_rep.json"
#     params:

#     message: "JSON: generate macs2 rep json"
#     log: _logfile
#     shell:
#         "cidc_chips/modules/scripts/json/json_macs2_rep.py  "
        
rule json_map:
    input:
        "analysis/align/{sample}/{sample}_mapping.txt"
    output:
        "analysis/json/{sample}/{sample}_map.json"
    message: "JSON: generate map json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_map.py -i {input} -o {output} "

rule json_meta:
    input:
        lambda wildcards: "analysis/ceas/%s/%s_summary.txt" % (getRunAndRep(wildcards),getRunAndRep(wildcards))
    output:
        "analysis/json/{sample}/{sample}_meta.json"
    message: "JSON: generate meta json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_meta.py -i {input} -o {output} -I {wildcards.sample} "

rule json_pbc:
    input:
        "analysis/frips/{sample}/{sample}_pbc.txt"
    output:
        "analysis/json/{sample}/{sample}_pbc.json"
    message: "JSON: generate pbc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_pbc.py -i {input} -o {output}"

# rule json_rep:
#     input:
#         rep_target
#     output:
#         "analysis/json/{sample}/{sample}_rep.json"
#     message: "JSON: generate rep json"
#     log: _logfile
#     shell:
#         "cidc_chips/modules/scripts/json/json_rep.py -c {input.cor} -o {output} -O {input.overlap} -i {run}" 


