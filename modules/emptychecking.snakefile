import os

def checking_targets(wildcards):
    ls = []
    ls.append(output_path + "/logs/empty_file_list.txt")
    return ls


def emptyCheckingInput(wildcards):
    # # Same as all_targets. Should be modified if all_targets modified.
    # # however DO NOT have checking_targets
    ls = all_targets(wildcards)
    for t in checking_targets(wildcards):
    	#print(t)
        if t in ls:
            ls.remove(t)
    return ls


rule emptyChecking_all:
    input:
        checking_targets

rule emptyChecking:
    input:
        emptyCheckingInput
    params:
        files = lambda wildards, input: " -f ".join(input),
    output:
        file=output_path + "/logs/empty_file_list.txt"
    message:
        "EMPTYCHECKING: checking whether any files are empty"
    shell:
        src_path + "/modules/scripts/emptycheck.py -f {params.files} -o {output}"













