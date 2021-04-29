# Author: Len Taing, Clara Cousins, Gali Bai
# Last modified: 01/11/2021
#MODULE: Chips report module

# Import packages
import yaml
from yaml import dump as yaml_dump

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
    
    #PEAK LEVEL QUALITY
    #ls.append(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.csv")
    ls.append(output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/02_details.yaml")
    
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

rule report_mapped_reads_bar:
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

#GALI TODO: write a helper script to calculate pbc using pbc.csv 
rule report_read_level_pbc:
    """Plot PBC"""
    input:
        output_path + "/frips/pbc.csv",
    output:
        csv=output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The PCR bottleneck coefficient (PBC) refers to the number of locations with exactly one uniquely mapped read divided by the number of unique genomic locations.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'PBC score','value':'PBC score'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""

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


#GALI TODO: write a helper script to combine the myco percentages
rule report_read_level_contamination_plot:
    """Plot contamination"""
    input:
        output_path + "/contam/contamination.csv"
    output:
        csv=output_path + "/report/Reads_Level_Quality/05_contamination_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/05_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The reported values for each species represent the percent of 100,000 reads that map to the reference genome of that species.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of 100,000 reads','value':'Percentage of 100,000 reads'}}}),
    group: "cohort_report"
    shell:
        #GALI TODO: contamination2.csv -> contamination2_parsed.csv??
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""

#Gali TODO- write a helper script to generate plotly line graph
# rule report_read_level_fragment_plot:
#     """Make fragment plots"""
#     input:
#         expand(output_path + "/frag/{sample}/{sample}_frags.txt", sample = list(config["samples"].keys()))
#     output:
#         output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly"
#     params:
#         files = lambda wildcards,inputs: " -f ".join(input) #A way to join input files for the helper script
#     shell:
#         pass
########################### END Read Level Quality Section ####################
########################### Peak Level Quality Section ########################
    
#GALI todo- write the shell part of this rule
#NOTE: it's ok to have derived report files as input
# rule report_peak_level_summary:
#     """Copy the csv files for rendering table of peak data"""
#     input:
#         output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
#         output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly",
#         output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
#         output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly"
#     output:
#         output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.csv"
#     shell:
#         pass

rule report_peak_level_peaks_plot:
    """Render number of peaks"""
    input:
        output_path + "/peaks/peakStats.csv",
    output:
        csv=output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The total peaks called, the peaks with a > 10 fold change (10FC), and the peaks with a > 20 fold change (20FC) for each run are represented here.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Number of peaks','value':'Number of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""

#GALI todo
# rule report_peaks_level_frips:
#     """Render FRIP"""
#     input:
#         #GALI to do--pull in the right input files for this plot
#     output:
#         csv=output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
#         details=output_path + "/report/Peaks_Level_Quality/03_details.yaml",
#     params:
#         files = lambda wildcards,input: " -f ".join(input),
#         caption="""caption: 'The fraction of reads in peaks (FRIP) score is the fraction of 4 million subsampled reads that fall within a defined peak region.'""",
#         plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'FRIP score (% of reads)','value':'FRIP score (% of reads)'}}}),
#     group: "cohort_report"
#     shell:
# #        """echo "{params.caption}" >> {output.details} &&
# #        echo "{params.plot_options}" >> {output.details} &&
# ##        cp {input} {output.csv}"""

#GALI TODO
# rule report_peaks_level_annotations:
#     input:
#         #GALI TODO- figure out which input files are needed to generate the output (when using a helper script to help parse)
#     output:
#         csv=output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
#         details=output_path + "/report/Peaks_Level_Quality/04_details.yaml",
#     params:
#         files = lambda wildcards,input: " -f ".join(input),
#         caption="""caption: 'The proportions of peaks for each sample overlapping with the promoters, exons, introns, and intergenic regions are shown here.'""",
#         plot_options = yaml_dump({'plotly': {'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of peaks','value':'Percentage of peaks'}}}),
#     group: "cohort_report"
#     shell:
#         pass

#GALI TODO
# rule report_peaks_level_dnase
# ...
########################### END Peak Level Quality Section ####################
rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_htmlTargets
    params:
        jinja2_template="cidc_chips/report/index.sample.html",
        output_path = output_path + "/report",
        #sections_list=",".join(['Overview', "Reads_Level_Quality", "Peaks_Level_Quality", "Genome_Track_View","Downstream"]), #define sections order here
        sections_list=",".join(['Overview','Reads_Level_Quality', 'Peaks_Level_Quality']), #define sections order here
        title="CHIPs Report",
    output:
        output_path+ "/report/report.html"
    message:
        "REPORT: Generating example report"
    shell:
        """python cidc_chips/modules/scripts/report.py -d {params.output_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r cidc_chips/report/static {params.output_path}"""
