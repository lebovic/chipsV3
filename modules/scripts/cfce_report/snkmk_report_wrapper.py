#!/usr/bin/env python
"""
Len Taing 2024 (TGBTG)

Script to call report.py--allows us to call from the shell: rather than 
through a run: in cfce_report.snakefile (rule cfce_report)
"""

import os
import sys
import argparse

from string import Template
#from chipsV3.modules.scripts.cfce_report.report import report
from report import report as snkmk_report
from report import data_uri

#DEPENDENCIES
from tabulate import tabulate

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

def processRunInfo(run_info_file):
    """extracts the macs version and the fdr used first and second line"""
    f = open(run_info_file)
    ver = f.readline().strip()
    fdr = f.readline().strip()
    return (ver,fdr)

def genPeakSummitsTable(conservPlots,motifSummary, runReps):
    """generates rst formatted table that will contain the conservation plot
    and motif analysis for ALL runs
    hdr = array of header/column elms
    NOTE: motifSummary is also acting like a boolean b/c we don't
    have access to config and need to replace if 'motif' in config...
    """
    #parse MotifSummary
    motifs = parseMotifSummary(motifSummary) if motifSummary else {}
    #runs = sorted(_getRepInput("$runRep"))
    runs = sorted(runReps)
    #HEADER- PROCESS the different modules differently
    #if 'motif' in config and config['motif'] == 'mdseqpos':
    if motifSummary and 'mdseqpos' in motifSummary:
        hdr = ["Run", "Conservation","MotifID","MotifName","Logo","Zscore"]
    else:
        hdr = ["Run", "Conservation","Motif","Logo","Pval","LogPval"]
    #BUILD up the rest of the table
    rest = []
    for run,img in zip(runs,conservPlots):
        #HANDLE null values
        if img and (img != 'NA'):
            conserv = ".. image:: %s" % data_uri(img)
        else:
            conserv = "NA"
        #HANDLE null values--Also check that we're doing motif analysis
        #if 'motif' in config and motifs[run]['logo'] and (motifs[run]['logo'] != 'NA'):
        if motifSummary and motifs[run]['logo'] and (motifs[run]['logo'] != 'NA'):
            motif_logo = ".. image:: %s" % data_uri(motifs[run]['logo'])
        else:
            motif_logo = "NA"

        #PROCESS the different modules differently
        #if 'motif' in config:
        if motifSummary:
            #if config['motif'] == 'mdseqpos':
            if 'mdseqpos' in motifSummary:
                rest.append([run, conserv, motifs[run]['motifId'], motifs[run]['motifName'], motif_logo,  motifs[run]['zscore']])
            else:
                rest.append([run, conserv, motifs[run]['motifName'], motif_logo, motifs[run]['pval'],motifs[run]['logp']])
        else:
            #motif analysis was skipped
            rest.append([run, conserv, 'NA', 'NA', 'NA','NA'])

    ret = tabulate(rest, hdr, tablefmt="rst")
    return ret

def parseMotifSummary(motif_csv):
    """Given a motifSummary.csv file, parses this into a dictionary 
    {run: {motifId: , motifName: , logo: , zscore: }}
    Returns this dictionary
    """
    ret = {}
    f = open(motif_csv)
    hdr = f.readline().strip().split(",")
    for l in f:
        tmp = l.strip().split(",")
        #if config['motif'] == 'mdseqpos':
        if 'mdseqpos' in motif_csv:
            ret[tmp[0]] = {'motifId': tmp[1], 'motifName': tmp[2], 'logo': tmp[3], 'zscore': tmp[4]}
        else:
            ret[tmp[0]] = {'motifName':tmp[1], 'logo':tmp[2], 'pval': tmp[3], 'logp': tmp[4]}
    f.close()
    return ret

def sampleGCandContam_table(fastqc_stats, fastqc_gc_plots, contam_table):
    """generates rst formatted table that will contain the fastqc GC dist. plot
    for all samples
    hdr = array of header/column elms
    NOTE: we use the thumb nail image for the GC content plots
    """
    #READ in fastqc_stats
    f_stats = {}
    f = open(fastqc_stats)
    hdr = f.readline().strip().split(",")  #largely ignored
    for l in f:
        tmp = l.strip().split(",")
        #store GC content, 3 col, using sample names as keys
        f_stats[tmp[0]] = tmp[2]
    f.close()
    
    #READ in contam_panel
    contam = {}
    f = open(contam_table)
    hdr = f.readline().strip().split(',') #We'll use this!
    species = hdr[1:]
    for l in f:
        tmp = l.strip().split(",")
        #build dict, use sample name as key
        contam[tmp[0]] = zip(species, tmp[1:]) #dict of tupes, (species, %)
    f.close()

    #PUT GC plots into a dictionary--keys are sample names
    plots = {}
    for p in fastqc_gc_plots:
        #NOTE the path structure is analysis/fastqc/{sample}/png_filename
        #we take second to last
        tmp = p.split("/")[-2]
        plots[tmp] = str(p)
    
    #build output
    ret=[]
    samples = sorted(list(f_stats.keys()))
    hdr = ["Sample", "GC median", "GC distribution"]
    hdr.extend(species)
    rest=[]

    for sample in samples:
        #HANDLE null values
        if plots[sample] and (plots[sample] != 'NA'):
            gc_plot = ".. image:: %s" % data_uri(plots[sample])
        else:
            gc_plot = "NA"
        #get rest of values and compose row
        contam_values = [v for (s, v) in contam[sample]]
        row = [sample, f_stats[sample], gc_plot]
        row.extend(contam_values)

        rest.append(row)
    ret = tabulate(rest, hdr, tablefmt="rst")
    return ret

def main():
    #usage = "USAGE: %prog --cfce_logo {filepath} --run_info {filepath} --map_stat {filepath} --pbc_stat {filepath} --peakFoldChange_png {filepath} --conservPlots {filepath} --samples_summary {filepath} --runs_summary {filepath} --contam_panel {filepath} --fastqc_stats {filepath} -p [fastqc_gc_plots 1] -p [fastqc_gc_plots 2] ... -[fastqc_gc_plots n] --motif {filepath}"
    parser = argparse.ArgumentParser(description="Call the legacy snakemake report generator")
    parser.add_argument("--cfce_logo", help="path to logo", required=True)
    parser.add_argument("--run_info", help="path to peaks/run_info.txt", required=True)
    parser.add_argument("--map_stat", help="path to mapping.png", required=True)
    parser.add_argument("--pbc_stat", help="path to pbc.png", required=True)
    parser.add_argument("--peakFoldChange_png", help="path to peakFoldChange.png", required=True)
    parser.add_argument("-c", "--conservPlots", action="append", help="path to runRep_conserv_thumb.png", required=True)
    parser.add_argument("--samples_summary", help="path to sequencingStatsSummary.csv", required=True)
    parser.add_argument("--runs_summary", help="path to peaksSummary.csv", required=True)
    parser.add_argument("--contam_panel", help="path to contamination_filtered.csv", required=True)
    parser.add_argument("--fastqc_stats", help="path to fastqc.csv", required=True)
    parser.add_argument("-f", "--fastq_gc_plots", action="append", help="paths to {sample}_perSeqGC_thumb.png", required=True)
    parser.add_argument("--motif", help="path to motifSummary.csv", default=None)
    #add argument for runs, src_path, output_path
    parser.add_argument("-r", "--runReps", action="append", help="chips run reps", required=True)  
    parser.add_argument("--src_path", help="chips src path", required=True)  
    parser.add_argument("--output_path", help="chips output path", required=True)  
    parser.add_argument("-o", "--output", help="output html file", required=True)

    args = parser.parse_args()

    output_path= args.output_path
    src_path= args.src_path
    runReps = args.runReps
    
    _cfce_report_log = output_path + "/logs/cfce_report.log"
    _CFCE_LOGO = src_path + "/static/cfce_report/CFCE_Logo_Final.jpg"
    _ReportCSS = src_path + "/static/cfce_report/report.css"
    _ReportTemplate = Template(open(src_path + "/static/cfce_report/cfce_report_template.txt").read())

    (macsVersion, fdr) = processRunInfo(args.run_info)
    #this local variable binds {samplesSummaryTable} in report template
    samplesSummaryTable = csvToSimpleTable(args.samples_summary)
    runsSummaryTable = csvToSimpleTable(args.runs_summary)
    #conservPlots = sorted(_getRepInput(output_path + "/conserv/$runRep/$runRep_conserv_thumb.png"))
    peakSummitsTable = genPeakSummitsTable(args.conservPlots, args.motif, runReps)
    #fastqc_gc_plots = [output_path + "/fastqc/%s/%s_perSeqGC_thumb.png" % (s,s) for s in config['samples']]
    sampleGCandContam = sampleGCandContam_table(args.fastqc_stats, args.fastq_gc_plots, args.contam_panel)

    #convert the args to dictionary, dropping any Nones
    #ALSO drop output, runs, src_path, output_path param b/c we don't want to embed that
    drops = ['output', 'output_path', 'src_path', 'runReps']
    args_dict = dict(filter(lambda x: x[1] and x[0] not in drops, args._get_kwargs()))
    #print(args_dict)
    tmp = _ReportTemplate.substitute(cfce_logo=data_uri(args.cfce_logo),map_stat=data_uri(args.map_stat),pbc_stat=data_uri(args.pbc_stat),peakFoldChange_png=data_uri(args.peakFoldChange_png),peakSummitsTable=peakSummitsTable)
    snkmk_report(tmp, args.output, stylesheet=_ReportCSS, metadata="Len Taing", **args_dict)

if __name__=='__main__':
    main()


