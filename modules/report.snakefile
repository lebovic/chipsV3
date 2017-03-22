#REPORT module
import sys
from string import Template
from snakemake.utils import report
from snakemake.report import data_uri

#DEPENDENCIES
from tabulate import tabulate

_ReportTemplate = Template(open("chips/static/chips_report.txt").read())
_logfile = "analysis/logs/report.log"

def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = ['analysis/report/samplesSummary.csv']
    ls.append('analysis/report/samplesSummary.csv')
    ls.append('analysis/report/runsSummary.csv')
    ls.append('analysis/report/report.html')
    return ls

def csvToSimpleTable(csv_file):
    """function to translate a .csv file into a reStructuredText simple table
    ref: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#simple-tables
    """
    #read in file
    f = open(csv_file)
    stats = [l.strip().split(",") for l in f]
    f.close()

    hdr = stats[0]
    rest = stats[1:]
    #RELY on the tabulate pkg
    ret = tabulate(rest, hdr, tablefmt="rst")
    return ret

rule report_all:
    input:
        report_targets

rule report:
    input:
        map_stat="analysis/align/mapping.png",
        pbc_stat="analysis/align/pbc.png",
	nonChrM_stat="analysis/frips/nonChrM_stats.png",
	samples_summary="analysis/report/samplesSummary.csv",
    output: html="analysis/report/report.html"
    run:
        samplesSummaryTable = csvToSimpleTable(input.samples_summary)
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

rule samples_summary_table:
    input:
        fastqc = "analysis/fastqc/fastqc.csv",
        mapping = "analysis/align/mapping.csv",
        pbc = "analysis/frips/pbc.csv",
        frag = "analysis/frag/fragSizes.csv",
        bam = "analysis/ceas/samples/bamRegionStats.csv"
    output:
        "analysis/report/samplesSummary.csv"
    log: _logfile
    shell:
        "chips/modules/scripts/get_sampleSummary.py -f {input.fastqc} -m {input.mapping} -p {input.pbc} -r {input.frag} -b {input.bam} > {output} 2>>{log}"

rule runs_summary_table:
    input:
        peaks = "analysis/peaks/peakStats.csv",
        frips = "analysis/frips/frips.csv",
    output:
        "analysis/report/runsSummary.csv"
    log: _logfile
    shell:
        "chips/modules/scripts/get_runsSummary.py -p {input.peaks} -f {input.frips} -o {output} 2>>{log}"
