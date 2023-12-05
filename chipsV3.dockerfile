FROM debian:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/envs/chips/share/homer/bin:$PATH
ENV PATH /opt/weblogo:$PATH

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git bc \
    libarchive-dev

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN apt-get install -y curl grep sed dpkg && \
    curl -L "https://github.com/krallin/tini/releases/download/v0.19.0/tini_0.19.0.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'

### GENERAL ###
#LEN NOTE: since we were getting
#"SyntaxError: future feature annotations is not defined" errors
#with python3.6.15 and markdown 3.5.1, I tried upgrading to python>=3.7.0
#but could not get the python to switch off of 3.6.15 (probably b/c snakemake)
#so I downgraded markdown to 3.4.4 and it works; also added dataclasses pkg
# DON'T change otherwise bug:
# docutils==0.16

RUN mamba create -n chipsV3 -c conda-forge -c bioconda snakemake=5.4.5 samtools>=1.10 python>=3.7.0 r-base=3.5.1 numpy>=1.19.5 pandas>=1.1.5 fastp=0.20.1 bwa=0.7.15 bowtie2=2.3.4.1 chromap>=0.2.5 sambamba=0.6.6 picard=2.18.4 bedtools=2.27.1 pybigwig>=0.3.17 seqtk>=1.3 fastqc>=0.11.5 r-ggplot2>=2.2.0 r-reshape2>=1.4.2 ruamel.yaml>=0.17.16 perl>=5.32.1 homer>=4.7 cython>=0.29.24 jinja2>=3.0.3 weblogo=2.8.2=3 bioconductor-seqlogo>=1.50.0 ghostscript ucsc-bedgraphtobigwig>=445 ucsc-bedsort>=445 macs2>=2.2.6 bioconductor-qdnaseq>=1.18.0 tabulate>=0.8.10 seaborn>=0.11.2 zip>=3.0  r-r.utils>=2.9.2 matplotlib>=3.1.1 intervaltree>=3.1.0 future>=0.18.2 hicmatrix>=16 pysam>=0.16.0 pytest>=6.2.5 pip>=21.3.1 gffutils>=0.12 pybedtools>=0.9.0 tqdm>=4.65.0 multiqc=1.9 plotly_express>=0.4.1 markdown=3.4.4 pyyaml>=5.4.1 pygenometracks>=3.6 docutils==0.16  dataclasses

#RUN apt-get install -y locales

#RUN locale-gen en_US.UTF-8

#COPY chipsDocker_env.sh /home
