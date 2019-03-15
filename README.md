# CPTAC3 RNA-seq splicing/transcript pipeline


## Setup
Create all the necessary packages inside a conda environment:

    conda create -n cptac3_splice \
        python=3.7 \
        ipython \
        snakemake-minimal \
        star=2.7 \
        samtools=1.9 htslib=1.9

Prepare a case list. For example, `case.list` is the default case list created by:

    # First 10 GBM samples in the Y2.b1 batch
    rg 'Y2.b1' /home/mwyczalk_test/Projects/CPTAC3/CPTAC3.catalog/CPTAC3.cases.dat \
        | rg 'GBM' \
        | cut -f1 \
        | head -n 10 > case.list

## Run the pipeline
Link all the required FASTQs:

    snakemake link_cptac_gdc_rna_fastqs

Run STAR alignment on all samples:

    snakemake star_align_all_samples
