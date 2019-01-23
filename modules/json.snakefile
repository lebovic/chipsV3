#MODULE: json- generating the analysis/json folder
_logfile="analysis/logs/json.log"

def json_targets(wildcards):
    ls = []
    ls.append("analysis/json")
    for sample in config["samples"]:
        ls.append("analysis/json/%s/%s_conserv.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_contam.json" % (sample, sample))
        if config["DHS"]:
            ls.append("analysis/json/%s/%s_dhs.json" % (sample, sample))
        if config["velcro_regions"]:
            ls.append("analysis/json/%s/%s_velcro.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_enrich_meta.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_fastqc.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_frag.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_frip.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_map.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_macs2_rep.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_macs2.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_meta.json" % (sample, sample))
        ls.append("analysis/json/%s/%s_rep.json" % (sample, sample))
    
    return ls

def comtamination_target(wildcards):
	return

def map_target(wildcards):
	return

def rep_target(wildcards):
	return

rule json_all:
    input:
        json_targets

rule json_conservation:
    input:
        "analysis/conserv/{run}.{rep}/{run}.{rep}_conserv.txt"
    output:
        "analysis/json/{sample}/{sample}_conserv.json"
    params:
        basics = "-b %s " % config["basics"] if config["basics"] else ""
        factor = "-f %s " % config["factor"] if config["factor"] else ""
        TF = "-T %s " % config["TF"] if config["TF"] else ""
    message: "JSON: generate conservation json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_conserv.py -i {input} -o {output} {params.basics} {params.factor} {params.TF} -I {run}"

rule json_comtamination:
    input:
        comtamination_target
    output:
        "analysis/json/{sample}/{sample}_contam.json"
    params:
        samples:
        id:
        species:
    message: "JSON: generate comtamination json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_comtamination.py  "

rule json_dhs:
    input:
        "analysis/ceas/{sample}/{sample}_DHS_summary.dhs"
    output:
        "analysis/json/{sample}/{sample}_dhs.json"
    message: "JSON: generate DHS json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_dhs.py -i {input} -o {output} -I {run}"

rule json_velcro:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_velcro.json"
    params:

    message: "JSON: generate velcro json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_velcro.py  "

rule json_enrich_meta:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_enrich_meta.json"
    params:

    message: "JSON: generate meta enrichment json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_enrich_meta.py  "

rule json_fastqc:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_fastqc.json"
    params:
        ids = "-s %s" % config[] if config[] else ""
    message: "JSON: generate fastqc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_fastqc.py -i {input} -o {output} {params.ids} -I {run}"

rule json_frag:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_frag.json"
    params:
        
    message: "JSON: generate frag json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frag.py  "
        
rule json_frip:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_frip.json"
    params:
        samples = "-s %s" % config[] if config[] else ""
    message: "JSON: generate frip json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frip.py -i {input} -o {output} {params.samples}"
        
rule json_macs2:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_macs2.json"
    message: "JSON: generate macs2 json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_macs2.py -i {input} -o {output} -I {run}"
        
rule json_macs2_rep:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_macs2_rep.json"
    params:

    message: "JSON: generate macs2 rep json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_macs2_rep.py  "
        
rule json_map:
    input:
        map_target
    output:
        "analysis/json/{sample}/{sample}_map.json"
    params:

    message: "JSON: generate map json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_map.py  "
        
rule json_meta:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_meta.json"
    message: "JSON: generate meta json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_meta.py -i {input} -o {output} -I {run} "
        
rule json_pbc:
    input:
        ""
    output:
        "analysis/json/{sample}/{sample}_pbc.json"
    params:
        samples = "-s %s" % config[] if config[] else ""
    message: "JSON: generate pbc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_pbc.py -i {input} -o {output} {parmas.samples} "
        
rule json_rep:
    input:
        rep_target
    output:
        "analysis/json/{sample}/{sample}_rep.json"
    message: "JSON: generate rep json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_rep.py -c {input.cor} -o {output} -O {input.overlap} -i {run}" 
        



