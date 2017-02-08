#REPORT module
from string import Template
from snakemake.utils import report
from snakemake.report import data_uri

_ReportTemplate = Template(open("chips/static/chips_report.txt").read())

rule report:
    input:
        map_stat="analysis/align/mapping.png",
        pbc_stat="analysis/align/pbc.png",
	nonChrM_stat="analysis/frips/nonChrM_stats.png"
    output: html="report.html"
    run:
        tmp = _ReportTemplate.substitute(map_stat=data_uri(input.map_stat),pbc_stat=data_uri(input.pbc_stat), nonChrM_stat=data_uri(input.nonChrM_stat))
        #report(_ReportTemplate, output.html, metadata="Len Taing", **input)
        report(tmp, output.html, metadata="Len Taing", **input)

rule plot_map_stat:
    input:
        "analysis/align/mapping.csv"
    output:
        "analysis/align/mapping.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/map_stats.R {input} {output}"

rule plot_pbc_stat:
    input:
        #"analysis/align/pbc.csv"
        "analysis/frips/pbc.csv"
    output:
        "analysis/align/pbc.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/plot_pbc.R {input} {output}"

rule plot_nonChrM_stats:
    input:
        "analysis/frips/nonChrM_stats.csv"
    output:
        "analysis/frips/nonChrM_stats.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/plot_nonChrM.R {input} {output}"
