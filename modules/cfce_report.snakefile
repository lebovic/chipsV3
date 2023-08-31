# Author: Len Taing (TGBTG) 2023
# Last modified: 2023-08-30
#MODULE: Legacy CFCE report module

#print(sys.path)
from string import Template
#from cidc_chips.modules.scripts.cfce_report.report import report
from scripts.cfce_report.report import report as snkmk_report
from scripts.cfce_report.report import data_uri

#DEPENDENCIES
from tabulate import tabulate

_CFCE_LOGO = "cidc_chips/static/cfce_report/CFCE_Logo_Final.jpg"
_ReportCSS = "cidc_chips/static/cfce_report/report.css"
_ReportTemplate = Template(open("cidc_chips/static/cfce_report/cfce_report_template.txt").read())

def cfce_report_targets(wildcards):
    ls = []
    ls.append(output_path + "/cfce_report/chips_report.html")
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

rule cfce_report_all:
    input:
        cfce_report_targets

def getReportInputs(wildcards):
    """Input function created just so we can switch-off motif analysis"""
    ret = {'cfce_logo': _CFCE_LOGO,
           'samples_summary': output_path + "/report/Reads_Level_Quality/01_read_level_summary_table.dt"}
    return ret

rule report:
    input:
        unpack(getReportInputs)
    output: html= output_path + "/cfce_report/chips_report.html"
    run:
        #this local variable binds {samplesSummaryTable} in report template
        samplesSummaryTable = csvToSimpleTable(input.samples_summary)
        tmp = _ReportTemplate.substitute(cfce_logo=data_uri(input.cfce_logo))
        snkmk_report(tmp, output.html, stylesheet=_ReportCSS, metadata="Len Taing", **input)
