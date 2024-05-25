#!/usr/bin/bash

#backup config.txt and data in $CONDA_PREFIX/share/homer/
mv $CONDA_PREFIX/share/homer/config.txt $CONDA_PREFIX/share/homer/config.txt.bak
mv $CONDA_PREFIX/share/homer/data $CONDA_PREFIX/share/homer/data.bak

#link in the homer files in ref_files
#NOTE: the conda dir is 3 levels from PROJECT, e.g. PROJECT/.snakemake/conda/{CONDA_PREFIX}

ln -s $CONDA_PREFIX/../../../ref_files/homer/config.txt $CONDA_PREFIX/share/homer/config.txt
ln -s $CONDA_PREFIX/../../../ref_files/homer/data $CONDA_PREFIX/share/homer/data
