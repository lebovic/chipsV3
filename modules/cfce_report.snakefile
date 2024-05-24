# Author: Len Taing (TGBTG) 2023
# Last modified: 2023-08-30
#MODULE: Legacy CFCE report module

#print(sys.path)
from string import Template
#from chipsV3.modules.scripts.cfce_report.report import report
from scripts.cfce_report.report import report as snkmk_report
from scripts.cfce_report.report import data_uri

#DEPENDENCIES
from tabulate import tabulate

_cfce_report_log = output_path + "/logs/cfce_report.log"
_CFCE_LOGO = src_path + "/static/cfce_report/CFCE_Logo_Final.jpg"
_ReportCSS = src_path + "/static/cfce_report/report.css"
_ReportTemplate = Template(open(src_path + "/static/cfce_report/cfce_report_template.txt").read())

def cfce_report_targets(wildcards):
    ls = []
    ls.append(output_path + "/cfce_report/cfce_chips_report.html")
    return ls

# moved to scripts/cfce_report/snkmk_report_wrapper.py
# def csvToSimpleTable(csv_file):
#     """function to translate a .csv file into a reStructuredText simple table
#     ref: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#simple-tables
#     """
#     #read in file
#     f = open(csv_file)
#     stats = [l.strip().split(",") for l in f]
#     f.close()

#     hdr = stats[0]
#     rest = stats[1:]
#     #RELY on the tabulate pkg
#     ret = tabulate(rest, hdr, tablefmt="rst")
#     return ret

# def processRunInfo(run_info_file):
#     """extracts the macs version and the fdr used first and second line"""
#     f = open(run_info_file)
#     ver = f.readline().strip()
#     fdr = f.readline().strip()
#     return (ver,fdr)

# def genPeakSummitsTable(conservPlots,motifSummary):
#     """generates rst formatted table that will contain the conservation plot
#     and motif analysis for ALL runs
#     hdr = array of header/column elms
#     """
#     #parse MotifSummary
#     motifs = parseMotifSummary(motifSummary) if motifSummary else {}
#     runs = sorted(_getRepInput("$runRep"))
#     #HEADER- PROCESS the different modules differently
#     if 'motif' in config and config['motif'] == 'mdseqpos':
#         hdr = ["Run", "Conservation","MotifID","MotifName","Logo","Zscore"]
#     else:
#         hdr = ["Run", "Conservation","Motif","Logo","Pval","LogPval"]
#     #BUILD up the rest of the table
#     rest = []
#     for run,img in zip(runs,conservPlots):
#         #HANDLE null values
#         if img and (img != 'NA'):
#             conserv = ".. image:: %s" % data_uri(img)
#         else:
#             conserv = "NA"
#         #HANDLE null values--Also check that we're doing motif analysis
#         if 'motif' in config and motifs[run]['logo'] and (motifs[run]['logo'] != 'NA'):
#             motif_logo = ".. image:: %s" % data_uri(motifs[run]['logo'])
#         else:
#             motif_logo = "NA"

#         #PROCESS the different modules differently
#         if 'motif' in config:
#             if config['motif'] == 'mdseqpos':
#                 rest.append([run, conserv, motifs[run]['motifId'], motifs[run]['motifName'], motif_logo,  motifs[run]['zscore']])
#             else:
#                 rest.append([run, conserv, motifs[run]['motifName'], motif_logo, motifs[run]['pval'],motifs[run]['logp']])
#         else:
#             #motif analysis was skipped
#             rest.append([run, conserv, 'NA', 'NA', 'NA','NA'])

#     ret = tabulate(rest, hdr, tablefmt="rst")
#     return ret

# def parseMotifSummary(motif_csv):
#     """Given a motifSummary.csv file, parses this into a dictionary 
#     {run: {motifId: , motifName: , logo: , zscore: }}
#     Returns this dictionary
#     """
#     ret = {}
#     f = open(motif_csv)
#     hdr = f.readline().strip().split(",")
#     for l in f:
#         tmp = l.strip().split(",")
#         if config['motif'] == 'mdseqpos':
#             ret[tmp[0]] = {'motifId': tmp[1], 'motifName': tmp[2], 'logo': tmp[3], 'zscore': tmp[4]}
#         else:
#             ret[tmp[0]] = {'motifName':tmp[1], 'logo':tmp[2], 'pval': tmp[3], 'logp': tmp[4]}
#     f.close()
#     return ret

# def sampleGCandContam_table(fastqc_stats, fastqc_gc_plots, contam_table):
#     """generates rst formatted table that will contain the fastqc GC dist. plot
#     for all samples
#     hdr = array of header/column elms
#     NOTE: we use the thumb nail image for the GC content plots
#     """
#     #READ in fastqc_stats
#     f_stats = {}
#     f = open(fastqc_stats)
#     hdr = f.readline().strip().split(",")  #largely ignored
#     for l in f:
#         tmp = l.strip().split(",")
#         #store GC content, 3 col, using sample names as keys
#         f_stats[tmp[0]] = tmp[2]
#     f.close()
    
#     #READ in contam_panel
#     contam = {}
#     f = open(contam_table)
#     hdr = f.readline().strip().split(',') #We'll use this!
#     species = hdr[1:]
#     for l in f:
#         tmp = l.strip().split(",")
#         #build dict, use sample name as key
#         contam[tmp[0]] = zip(species, tmp[1:]) #dict of tupes, (species, %)
#     f.close()

#     #PUT GC plots into a dictionary--keys are sample names
#     plots = {}
#     for p in fastqc_gc_plots:
#         #NOTE the path structure is analysis/fastqc/{sample}/png_filename
#         #we take second to last
#         tmp = p.split("/")[-2]
#         plots[tmp] = str(p)
    
#     #build output
#     ret=[]
#     samples = sorted(list(f_stats.keys()))
#     hdr = ["Sample", "GC median", "GC distribution"]
#     hdr.extend(species)
#     rest=[]

#     for sample in samples:
#         #HANDLE null values
#         if plots[sample] and (plots[sample] != 'NA'):
#             gc_plot = ".. image:: %s" % data_uri(plots[sample])
#         else:
#             gc_plot = "NA"
#         #get rest of values and compose row
#         contam_values = [v for (s, v) in contam[sample]]
#         row = [sample, f_stats[sample], gc_plot]
#         row.extend(contam_values)

#         rest.append(row)
#     ret = tabulate(rest, hdr, tablefmt="rst")
#     return ret

def getReportInputs(wildcards):
    """Input function created just so we can switch-off motif analysis"""
    ret = {'cfce_logo': _CFCE_LOGO,
           'run_info':output_path + "/peaks/run_info.txt",
           'map_stat':output_path + "/cfce_report/mapping.png",
           'pbc_stat':output_path + "/cfce_report/pbc.png",
           'peakFoldChange_png':output_path + "/cfce_report/peakFoldChange.png",
           'conservPlots': sorted(_getRepInput(output_path + "/conserv/$runRep/$runRep_conserv_thumb.png")),
           'samples_summary': output_path + "/cfce_report/sequencingStatsSummary.csv",
           'runs_summary':output_path + "/cfce_report/peaksSummary.csv",
           'contam_panel':output_path + "/contam/contamination_filtered.csv",
           'fastqc_stats':output_path + "/fastqc/fastqc.csv",
           'fastqc_gc_plots': expand(output_path + "/fastqc/{sample}/{sample}_perSeqGC_thumb.png", sample=config["samples"])
           }
    #MOTIF handled outside so we can check if it's run or not
    if 'motif' in config:
        ret['motif'] = output_path + "/motif/motifSummary.csv"
    return ret

rule cfce_report_all:
    input:
        cfce_report_targets

rule cfce_report:
    input:
        unpack(getReportInputs)
    params:
        conservPlots=lambda wildcards, input: " -c ".join(input.conservPlots),
        fastqcPlots=lambda wildcards,input: " -f ".join(input.fastqc_gc_plots),
        motif = lambda wildcards,input: "--motif %s" % input.motif if 'motif' in input else "",
        output_path= "--output_path " + output_path,
        src_path= "--src_path " + src_path,
        runs = " -r ".join(_getRepInput("$runRep")),
    conda: "../envs/report/cfce_report.yaml"
    output: html= output_path + "/cfce_report/cfce_chips_report.html"
    shell:
        src_path + "/modules/scripts/cfce_report/snkmk_report_wrapper.py --cfce_logo {input.cfce_logo} --run_info {input.run_info} --map_stat {input.map_stat} --pbc_stat {input.pbc_stat} --peakFoldChange_png {input.peakFoldChange_png} -c {params.conservPlots} --samples_summary {input.samples_summary} --runs_summary {input.runs_summary} --contam_panel {input.contam_panel} --fastqc_stats {input.fastqc_stats} -f {params.fastqcPlots} {params.motif} {params.output_path} {params.src_path} -r {params.runs} -o {output}"
    # run:
    #     (macsVersion, fdr) = processRunInfo(input.run_info)
    #     #this local variable binds {samplesSummaryTable} in report template
    #     samplesSummaryTable = csvToSimpleTable(input.samples_summary)
    #     runsSummaryTable = csvToSimpleTable(input.runs_summary)
    #     conservPlots = sorted(_getRepInput(output_path + "/conserv/$runRep/$runRep_conserv_thumb.png"))
    #     if 'motif' in config:
    #         peakSummitsTable = genPeakSummitsTable(conservPlots, input.motif)
    #     else:
    #         peakSummitsTable = genPeakSummitsTable(conservPlots, None)
    #     #HACK b/c unpack is ruining the input.fastqc_gc_plots element--the list
    #     #becomes a singleton
    #     fastqc_gc_plots = [output_path + "/fastqc/%s/%s_perSeqGC_thumb.png" % (s,s) for s in config['samples']]
    #     sampleGCandContam = sampleGCandContam_table(input.fastqc_stats, fastqc_gc_plots, input.contam_panel)
    #     tmp = _ReportTemplate.substitute(cfce_logo=data_uri(input.cfce_logo),map_stat=data_uri(input.map_stat),pbc_stat=data_uri(input.pbc_stat),peakFoldChange_png=data_uri(input.peakFoldChange_png),peakSummitsTable=peakSummitsTable)
    #     snkmk_report(tmp, output.html, stylesheet=_ReportCSS, metadata="Len Taing", **input)

rule samples_summary_table:
    input:
        fastqc = output_path + "/fastqc/fastqc.csv",
        mapping = output_path + "/align/mapping.csv",
        pbc = output_path + "/frips/pbc.csv",
    output:
        output_path + "/cfce_report/sequencingStatsSummary.csv"
    log: _cfce_report_log
    shell:
        src_path + "/modules/scripts/get_sampleSummary.py -f {input.fastqc} -m {input.mapping} -p {input.pbc} > {output} 2>>{log}"

rule runs_summary_table:
    input:
        peaks = output_path + "/peaks/peakStats.csv",
        frips = output_path + "/frips/frips.csv",
        dhs = output_path + "/ceas/dhs.csv",
        meta = output_path + "/ceas/meta.csv",
    output:
        output_path + "/cfce_report/peaksSummary.csv"
    log: _cfce_report_log
    shell:
        src_path + "/modules/scripts/get_runsSummary.py -p {input.peaks} -f {input.frips} -d {input.dhs} -m {input.meta} -o {output} 2>>{log}"

######## PLOTS ######
rule plot_map_stat:
    input:
        output_path + "/align/mapping.csv"
    output:
        output_path + "/cfce_report/mapping.png"
    conda: "../envs/report/cfce_report.yaml"
    log: _cfce_report_log
    shell:
        "Rscript " + src_path + "/modules/scripts/report_mapStats.R {input} {output}"

rule plot_pbc_stat:
    input:
        output_path + "/frips/pbc.csv"
    output:
        output_path + "/cfce_report/pbc.png"
    conda: "../envs/report/cfce_report.yaml"
    log: _cfce_report_log
    shell:
        "Rscript " + src_path + "/modules/scripts/report_plotPBC.R {input} {output}"

rule plot_peakFoldChange:
    input: 
        output_path + "/peaks/peakStats.csv"
    output:
        output_path + "/cfce_report/peakFoldChange.png"
    conda: "../envs/report/cfce_report.yaml"
    log: _cfce_report_log
    shell:
        "Rscript " + src_path + "/modules/scripts/report_plotFoldChange.R {input} {output}"
