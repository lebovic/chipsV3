#!/usr/bin/env python
"""
Len Taing 2024 (TGBTG)

Script to call the stuff that was formerly in motif_homer rule's run section
"""

import os
import sys
import argparse
import subprocess

def _createEmptyMotif(motif_html):
    """When the _sorted_5k_summits.bed has too few peaks, or is empty,
    we still want to create an emtpy homerResult.html
    INPUT: output paths of these files
    """
    #CHECK for dir existence:
    _path = "/".join(motif_html.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    #Create an empty mdseqpos_index.html
    subprocess.call(['touch', motif_html])


def main():
    parser = argparse.ArgumentParser(description="Call the motif homer IF the bed file has more than _minPeaks, otherwise create empty outputs")
    parser.add_argument("-b", "--bed", help="path to bed file", required=True)
    parser.add_argument("-m", "--minPeaks", help="minimum number of peaks (default: 500)", default="500")
    parser.add_argument("-g", "--genome", help="genome", required=True)
    parser.add_argument("-r", "--results_path", help="output results path", required=True)
    parser.add_argument("-s", "--size", help="motif window size", default=600)
    parser.add_argument("-t", "--threads", help="number of threads", default=8)
    parser.add_argument("--output_html", help="html output path", required=True)
    parser.add_argument("--output_txt", help="txt output path", required=True)
    parser.add_argument("--output_logo", help="logo output path", required=True)

    args = parser.parse_args()

    wc = str(subprocess.check_output(['wc', '-l', args.bed]))
    #this returns a byte string--we need to eliminate the b' ... '
    #and then convert the first elm to int
    wc = int(wc[2:-1].split()[0])

    if wc >= int(args.minPeaks):
        try:
            #PASS- run motif scan
            cmd = ["findMotifsGenome.pl", args.bed, args.genome, args.results_path, "-size", args.size, "-p", args.threads, "-mask", "-seqlogo", "-preparsedDir", args.results_path]
            subprocess.call(cmd)
            #shell("findMotifsGenome.pl {input} {params.genome} {params.results} -size {params.size} -p {threads} -mask -seqlogo -preparsedDir {params.results} >>{log} 2>&1")
        except:
            print("homer package not installed")
    else:
        #FAIL - create empty outputs
        _createEmptyMotif(args.output_html)
        _createEmptyMotif(args.output_txt)
        _createEmptyMotif(args.output_logo)

if __name__=='__main__':
    main()


