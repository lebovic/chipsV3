_fasp_threads=2

def trim_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	# fastp output json and html summarization
    	ls.append(output_path + "/trim_adaptor/%s/%s_fastp.json" % (sample,sample))
    	ls.append(output_path + "/trim_adaptor/%s/%s_fastp.html" % (sample,sample))
    	fastq_num = len(config["samples"][sample]) 
    	if fastq_num > 1: #paired-end
    		ls.append(expand(output_path + "/trim_adaptor/%s/%s_{mate}.trimmed.fq" % (sample,sample), mate = range(fastq_num)))
    	else: # single-end
    		ls.append(output_path + "/trim_adaptor/%s/%s_1.trimmed.fq" % (sample,sample))
    return ls

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

## it is not possible to write an output function to determine the output file number
## based on the wildcards. use touch to trick snakemake instead.
## see https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
## and https://stackoverflow.com/questions/56861913/split-bam-by-clusters-and-then-merge-bam-by-cluster-using-checkpoint

FASTP_VERSION = subprocess.check_output("fastp -v", shell=True)

rule trim_fasp:
	input: getFastq
	output: 
		fastqs = expand(output_path + "/trim_adaptor/{{sample}}/{{sample}}_{mate}.trimmed.fq", mate = range(2)),
		json = output_path + "/trim_adaptor/{sample}/{sample}_fastp.json",
		html = output_path + "/trim_adaptor/{sample}/{sample}_fastp.html" 
	threads: _fasp_threads
	message: "trimming adaptors for {input} using fastp"
	log: output_path + "/logs/trim_adaptor/{sample}.log"
	conda: "../envs/trimming/trim_fastp.yaml"
	run:
		if len(input) > 1:
			shell("fastp --thread {threads} --detect_adapter_for_pe --in1 {input[0]} --in2 {input[1]} --out1 {output.fastqs[0]} --out2 {output.fastqs[1]} -h {output.html} -j {output.json} > {log} 2>&1 ")
		else:
			shell("fastp --thread {threads} --in1 {input[0]} --out1 {output[0]} -h {output.html} -j {output.json} > {log} 2>&1")
			shell("touch {output[1]}")
	
