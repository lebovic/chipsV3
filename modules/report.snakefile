#REPORT module
from string import Template
from snakemake.utils import report
from snakemake.report import data_uri

_ReportTemplate = Template(open("chips/static/chips_report.txt").read())

rule report:
    input:
        map_stat="analysis/align/mapping.png"
    output: html="report.html"
    run:
        tmp = _ReportTemplate.substitute(map_stat=data_uri(input.map_stat))
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
