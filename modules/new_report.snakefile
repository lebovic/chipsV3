# Author: Len Taing, Clara Cousins, Gali Bai
# Last modified: 01/11/2021
#MODULE: Chips report module

# Import packages
import yaml
from yaml import dump as yaml_dump
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import re
import glob

def report_htmlTargets(wildcards):
    ls = []
    ls.append(output_path + "/report/Overview/01_chips_workflow.png")
    ls.append(output_path + "/report/Overview/01_details.yaml")
    ls.append(output_path + "/report/Overview/02_select_software_versions.tsv")
    ls.append(output_path + "/report/Overview/02_details.yaml")
    ls.append(output_path + "/report/Overview/03_assembly.csv")

    #READ LEVEL QUALITY
    ls.append(output_path + "/report/Reads_Level_Quality/01_read_level_summary_table.dt")
    ls.append(output_path + "/report/Reads_Level_Quality/01_details.yaml")
    ls.append(output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/02_details.yaml")
    ls.append(output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/03_details.yaml")
    ls.append(output_path + "/report/Reads_Level_Quality/04_contamination_table.dt")
    ls.append(output_path + "/report/Reads_Level_Quality/04_details.yaml")
    ls.append(output_path + "/report/Reads_Level_Quality/05_contamination_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/05_details.yaml")
    ls.append(output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/06_details.yaml")

    #PEAK LEVEL QUALITY
    ls.append(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.dt")
    ls.append(output_path + "/report/Peaks_Level_Quality/01_details.yaml")
    ls.append(output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/02_details.yaml")
    ls.append(output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/03_details.yaml")
    ls.append(output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/04_details.yaml")
    ls.append(output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/05_details.yaml")

    #Genome Track View
    for list_num, gene in enumerate(config["genes_to_plot"].strip().split()):
        ls.append((output_path + "/report/Genome_Track_View/{num}_genome_track_for_{track}.png").format(num = list_num, track = gene))
    return ls

rule report_all:
    input:
        output_path+ "/report/report.html"

########################### OVERVIEW Section ##################################
rule report_overview_workflow:
    input:
        png="cidc_chips/report/chips_workflow.png",
        yml="cidc_chips/report/intro_details.yaml",
    output:
        png = output_path + "/report/Overview/01_chips_workflow.png",
        det = output_path + "/report/Overview/01_details.yaml",
    shell:
        "cp {input.png} {output.png} && cp {input.yml} {output.det}"

rule report_overview_software_versions:
    #inpuy: #NO Input
    output:
        tsv=output_path + "/report/Overview/02_select_software_versions.tsv",
        details=output_path + "/report/Overview/02_details.yaml",
    params:
        caption="""caption: 'The details of other software used in CHIPs are written to software_versions_all.txt in the directory where this report was generated.'"""
    shell:
        """echo "{params.caption}" >> {output.details} && """
        """cidc_chips/modules/scripts/report/overview/software_versions.py -o {output.tsv}"""

rule report_overview_assembly:
    #inpuy: #NO Input
    output:
        output_path + "/report/Overview/03_assembly.csv",
    params:
        assembly=config['assembly']
    shell:
        "echo {params.assembly} > {output}"
########################### END OVERVIEW Section ##############################

########################### Read Level Quality Section ########################
rule report_read_level_summary_table:
    input:
        mapping=output_path + "/align/mapping.csv",
        pbc=output_path + "/frips/pbc.csv",
    output:
        csv=output_path + "/report/Reads_Level_Quality/01_read_level_summary_table.dt",
        details=output_path + "/report/Reads_Level_Quality/01_details.yaml",
    params:
        caption="""caption: 'Abbreviations: M, million; PBC, PCR bottlneck coefficient.'"""
    shell:
        """echo "{params.caption}" >> {output.details} && """
        """cidc_chips/modules/scripts/report/read_level_quality/read_level_summary.py -m {input.mapping} -p {input.pbc} -o {output.csv}"""

rule report_read_level_mapped_reads:
    input:
        output_path + "/align/mapping.csv",
    output:
        csv=output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/02_details.yaml",
    params:
        caption="""caption: 'Mapped reads refer to the number of reads successfully mapping to the genome, while uniquely mapped reads are the subset of mapped reads mapping only to one genomic location.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Number of reads','value':'Number of reads'}}}),
    shell:
        """echo "{params.caption}" >> {output.details} && """
        """echo "{params.plot_options}" >> {output.details} &&"""
        """cp {input} {output.csv}"""

rule report_read_level_pcr_bottleneck_coefficient:
    """Plot PBC"""
    input:
        output_path + "/frips/pbc.csv",
    output:
        csv=output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/03_details.yaml",
    params:
        caption="""caption: 'The PCR bottleneck coefficient (PBC) refers to the number of locations with exactly one uniquely mapped read divided by the number of unique genomic locations.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'PBC score','value':'PBC score'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cidc_chips/modules/scripts/report/read_level_quality/read_level_pbc.py -p {input} -o {output.csv}"""

rule report_read_level_contamination_tbl:
    input:
         output_path + "/contam/contamination.csv"
    output:
         csv=output_path + "/report/Reads_Level_Quality/04_contamination_table.dt",
         details=output_path + "/report/Reads_Level_Quality/04_details.yaml"
    params:
         caption="caption: 'Contamination percentages for all reference genomes are included here.' "
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cp {input} {output.csv}"""

rule report_read_level_contamination_plot:
    """Plot contamination"""
    input:
        output_path + "/contam/contamination.csv"
    output:
        csv=output_path + "/report/Reads_Level_Quality/05_contamination_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/05_details.yaml",
    params:
        #files=lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The reported values for each species represent the percent of 100,000 reads that map to the reference genome of that species.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of 100,000 reads','value':'Percentage of 100,000 reads'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cidc_chips/modules/scripts/report/read_level_quality/read_level_contam.py -c {input} -o {output.csv}"""

rule report_read_level_fragment_plot:
    """Make fragment plots"""
    input:
        expand(output_path + "/frag/{sample}/{sample}_frags.txt", sample = list(config["samples"].keys())),
    output:
        csv=output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly",
        details=output_path + "/report/Reads_Level_Quality/06_details.yaml",
    params:
        files=lambda wildcards, input: " -f ".join(input),
        caption="""caption: 'Fragment size distributions show paired-end fragments in each sample. The plotted value for each sample is the probability density in a 5 bp bin size normalized so the integral is 1.'""",
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cidc_chips/modules/scripts/report/read_level_quality/read_level_frag.py -f {params.files} -o {output.csv}"""

########################### END Read Level Quality Section ####################
########################### Peak Level Quality Section ########################

rule report_peak_level_summary:
    """Copy the csv files for rendering table of peak data"""
    input:
        map= output_path + "/align/mapping.csv",
        peak= output_path + "/peaks/peakStats.csv",
        frip= output_path + "/frips/frips.csv",
        ceas= output_path + "/ceas/meta.csv",
        dhs= output_path + "/ceas/dhs.csv",
    output:
        details= output_path + "/report/Peaks_Level_Quality/01_details.yaml",
        sum= output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.dt",
        frip= output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
        ceas= output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
        dhs= output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly",
    params:
        caption="caption: 'Abbreviations: 10FC, > 10 fold change; 20FC, > 20 fold change; FRiP, Fraction of reads in peaks; Prom, Promoter; Inter, Intergenic; DHS, DNAseI hypersensitivity sites' "
    shell:
        """
        echo "{params.caption}" >> {output.details} &&
        cidc_chips/modules/scripts/report/peak_level_quality/peak_level_summary.py -p {input.peak} -f {input.frip} -m {input.ceas} -d {input.dhs} -s {output.sum} -r {output.frip} -a {output.ceas} -o {output.dhs}"""

rule report_peak_level_peaks_plot:
    """Render number of peaks"""
    input:
        output_path + "/peaks/peakStats.csv",
    output:
        csv=output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/02_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The total peaks called, the peaks with a > 10 fold change (10FC), and the peaks with a > 20 fold change (20FC) for each run are represented here.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Number of peaks','value':'Number of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""


rule report_peaks_level_frips:
    """Render FRIP"""
    output:
        #csv=output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/03_details.yaml",
    params:
        caption="""caption: 'The fraction of reads in peaks (FRIP) score is the fraction of 4 million subsampled reads that fall within a defined peak region.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'FRIP score (% of reads)','value':'FRIP score (% of reads)'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details}"""

rule report_peaks_level_annotations:
    """Render peak annotation"""
    output:
        #csv=output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/04_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The proportions of peaks for each sample overlapping with the promoters, exons, introns, and intergenic regions are shown here.'""",
        plot_options = yaml_dump({'plotly': {'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of peaks','value':'Percentage of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details}"""

rule report_peaks_level_dhs:
    """Render peak dhs"""
    output:
        #csv=output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/05_details.yaml",
    params:
        caption="""caption: 'DNAse hypersensitive sites (DHS) may represent highly active regions of the genome. The data below represent the percentage of 4 million subsampled peaks that intersect with DHS peaks as defined by list of known DHS regions (specific to each species).'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of peaks','value':'Percentage of peaks'}}}),
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details}"""

########################### END Peak Level Quality Section ####################
########################### Genome Track View Section ########################
rule report_genome_track_make_bed:
    """Make new bed files for genome track"""
    input:
        config['geneBed'],
    output:
        extend= output_path + "/report/Genome_Track_View/extend.bed",
        tss= output_path + "/report/Genome_Track_View/tss.bed"
    params:
        up= config['upstream'],
        down= config['downstream'],
    shell:
        """cidc_chips/modules/scripts/report/genome_track_view/make_bed_file.py -i {input} -u {params.up} -d {params.down} -e {output.extend} -t {output.tss}"""

def genome_tracks_init_inputFn(wildcards):
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ls.append((output_path + "/peaks/{runRep}/{runRep}_treat_pileup.bw").format(runRep = runRep))
    return ls

rule report_genome_track_make_tracks:
    """Make genome track configuration file"""
    input:
        genome_tracks_init_inputFn
    output:
        output_path + "/report/Genome_Track_View/tracks_all_vlines.ini",
    params:
        track= temp(output_path + "/report/Genome_Track_View/tracks_all.ini"),
        extend= output_path + "/report/Genome_Track_View/extend.bed",
        tss= output_path + "/report/Genome_Track_View/tss.bed",
    shell:
        """make_tracks_file --trackFiles {input} -o {params.track} &&
        cidc_chips/modules/scripts/report/genome_track_view/make_track_file.py -i {params} -e {params.extend} -t {params.tss} -o {output}"""

png_list = []
for list_num, gene in enumerate(config["genes_to_plot"].strip().split()):
    png_list.append((output_path + "/report/Genome_Track_View/{num}_genome_track_for_{track}.png").format(num = list_num, track = gene))


rule report_genome_track_make_plot:
    """Make genome track plot"""
    input:
        ini= output_path + "/report/Genome_Track_View/tracks_all_vlines.ini",
        extend= output_path + "/report/Genome_Track_View/extend.bed",
    output:
        png_list
        #png= lambda wildcards: [" -o " + output_path + "/report/Genome_Track_View/%s_genome_track_for_%s.png" % (i,i) for i, i in enumerate(config["genes_to_plot"].strip().split())]
    params:
        genes= lambda wildcards: [" -g %s" % g for g in config["genes_to_plot"].strip().split()],
        png= lambda wildcards, output: " -o ". join(output),
    shell:
        """cidc_chips/modules/scripts/report/genome_track_view/make_track_png.py -i {input.ini} -e {input.extend} {params.genes} -o {params.png}"""

########################### END Genome Track View Section ####################
########################### Downstreat Section ################################
def report_downstream_conser_motif_inputFn(wildcards):
    #get conservation png files
    runRep_ls = []
    conserv_ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            runRep_ls.append(runRep)
            conserv_ls.append((output_path + "/conserv/{runRep}/{runRep}_conserv_thumb.png").format(runRep = runRep))

    #get homer ls files
    homer_ls = []
    if config['motif'] == 'homer':
        for run in config["runs"].keys():
            for rep in _reps[run]:
                runRep = "%s.%s" % (run, rep)
                #LEN: Please check this path!
                conserv_ls.append((output_path + "/motif/{runRep}/results/knownResults/known1.logo.pnt").format(runRep = runRep))
    tmp = {'conserv_logos': conserv_ls,
           'homer_logos': homer_ls,
           'runReps': runRep_ls}
    return tmp

rule report_downstream_conser_motif:
    input:
        unpack(report_downstream_conser_motif_inputFn)
    output:
        csv = output_path + "/report/Downstream/01_conservation_and_top_motifs.csv",
        det = output_path + "/report/Downstream/01_details.yaml"
    params:
        outpath="report/Downstream/",
        conserv_logos= lambda wildcards, input: " -c ".join(input.conserv_logos),
        motif_logos= lambda wildcards, input: " -m ".join(input.motif_logos),
        runReps = lambda wildcards, input: " -r ".join(input.runReps),
    shell:
        #TODO:need to look up name and log pval for each runRep's homer results
        "cidc_chips/modules/scripts/report/downstream/conserv_motif_table.py -r {params.runReps} -c {params.conserv_logos} -m {params.motif_logos}"
    
########################### END Downstreat Section ############################
rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_htmlTargets
    params:
        jinja2_template="cidc_chips/report/index.sample.html",
        output_path = output_path + "/report",
        #sections_list=",".join(['Overview', "Reads_Level_Quality", "Peaks_Level_Quality", "Genome_Track_View","Downstream"]), #define sections order here
        sections_list=",".join(['Overview','Reads_Level_Quality', 'Peaks_Level_Quality', 'Genome_Track_View']), #define sections order here
        title="CHIPs Report",
    output:
        output_path+ "/report/report.html"
    message:
        "REPORT: Generating example report"
    shell:
        """python cidc_chips/modules/scripts/report.py -d {params.output_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r cidc_chips/report/static {params.output_path}"""
