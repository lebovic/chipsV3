#MODULE: json- generating the analysis/json folder
_logfile="analysis/logs/json.log"

def json_targets(wildcards):
    ls = []
    ls.append("analysis/json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_conserv.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_contam.json")
    if :
        ls.append("analysis/json/{run}.{rep}/{run}.{rep}_dhs.json")
    if :
        ls.append("analysis/json/{run}.{rep}/{run}.{rep}_velcro.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_enrich_meta.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_fastqc.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_frag.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_frip.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_map.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_macs2_rep.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_macs2.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_meta.json")
    ls.append("analysis/json/{run}.{rep}/{run}.{rep}_rep.json")
    
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
        "analysis/json/{run}.{rep}/{run}.{rep}_conserv.json"
    params:
        basics = "-b %s " % config[]
        factor = "-f %s " % config[]
        TF = "-T %s " % config[]
    message: "JSON: generate conservation json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_conserv.py -i {input} -o {output} {params.basics} {params.factor} {params.TF} -I {run}"

rule json_comtamination:
    input:
        comtamination_target
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_contam.json"
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
        "analysis/ceas/{run}.{rep}/{run}.{rep}_DHS_summary.dhs"
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_dhs.json"
    message: "JSON: generate DHS json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_dhs.py -i {input} -o {output} -I {run}"

rule json_velcro:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_velcro.json"
    params:

    message: "JSON: generate velcro json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_velcro.py  "

rule json_enrich_meta:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_enrich_meta.json"
    params:

    message: "JSON: generate meta enrichment json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_enrich_meta.py  "

rule json_fastqc:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_fastqc.json"
    params:
        ids = "-s %s" % config[]
    message: "JSON: generate fastqc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_fastqc.py -i {input} -o {output} {params.ids} -I {run}"

rule json_frag:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_frag.json"
    params:
        
    message: "JSON: generate frag json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frag.py  "
        
rule json_frip:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_frip.json"
    params:
        samples = "-s %s" % config[]
    message: "JSON: generate frip json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_frip.py -i {input} -o {output} {params.samples}"
        
rule json_macs2:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_macs2.json"
    message: "JSON: generate macs2 json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_macs2.py -i {input} -o {output} -I {run}"
        
rule json_macs2_rep:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_macs2_rep.json"
    params:

    message: "JSON: generate macs2 rep json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_macs2_rep.py  "
        
rule json_map:
    input:
        map_target
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_map.json"
    params:

    message: "JSON: generate map json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_map.py  "
        
rule json_meta:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_meta.json"
    message: "JSON: generate meta json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_meta.py -i {input} -o {output} -I {run} "
        
rule json_pbc:
    input:
        ""
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_pbc.json"
    params:
        samples = "-s %s" % config[]
    message: "JSON: generate pbc json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_pbc.py -i {input} -o {output} {parmas.samples} "
        
rule json_rep:
    input:
        rep_target
    output:
        "analysis/json/{run}.{rep}/{run}.{rep}_rep.json"
    message: "JSON: generate rep json"
    log: _logfile
    shell:
        "cidc_chips/modules/scripts/json/json_rep.py -c {input.cor} -o {output} -O {input.overlap} -i {run}" 
        



