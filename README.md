# CPTAC3 RNA-seq splicing/transcript pipeline


## Setup
Create all the necessary packages inside a conda environment:

    conda create -n cptac3_splice \
        python=3.7 \
        ipython \
        snakemake-minimal


## Run the pipeline

    snakemake ...
