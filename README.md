# CPTAC3 RNA-seq splicing/transcript pipeline



## Setup
Create all the necessary packages inside a conda environment:

    conda create -n cptac3_splice \
        python=3.7 \
        ipython \
        snakemake-minimal \
        star=2.7 stringtie=1.3 \
        samtools=1.9 htslib=1.9 \
        mapsplice=2.2 \
        gffcompare=0.10

Prepare a case list. For example, `case.list` is the default case list created by:

    # First 10 GBM samples in the Y2.b1 batch
    rg 'Y2.b1' /home/mwyczalk_test/Projects/CPTAC3/CPTAC3.catalog/CPTAC3.cases.dat \
        | rg 'GBM' \
        | cut -f1 \
        | head -n 10 > case.list



## Run the pipeline
The pipeline defines two resources:
- `mem_mb`: Memory limit in MB
- `io_heavy`: Maximal number of concurrent IO heavy tasks

Specify the limits while running any snakemake job. For example, to use 100GB
of memory and 4 concurrent IO heavy tasks using 20 CPU cores,

    snakemake -j20 --resources mem_mb=100000 io_heavy=4 ...

List all the available commands by

    snakemake -l


### Sample setup
Link all the required FASTQs:

    snakemake link_cptac_gdc_rna_fastqs


### StringTie
Run STAR alignment on all samples:

    snakemake star_align_all_samples

Run StringTie novel transcript discovery of all samples:

    snakemake stringtie_merge_gtfs

Transcript GTFs will be available at:
- Per sample at `processed_data/stringtie/{sample}/transcripts.gtf`
- All samples merged at `processed_data/stringtie/merged.gtf`


### MapSplice2
MapSplice requires the reference chromosomes to be stored in individual FASTA
files, which can be done by the script at `scripts/split_genome_fa_per_chrom.sh`:

    bash scripts/split_genome_fa_per_chrom.sh

By default the per chromosome FASTAs are stored at
`/diskmnt/Projects/cptac_scratch/GRCh38.d1.vd1_per_chrom`.

MapSplice requires Bowtie1 index for alignment. Since MapSplice can only take
FASTA witout any other annotation in the header, we strip all the annotations
in the FASTA file first by

    sed -r 's/^(>[^ ]+) .*$/\1/' \
        /diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa \
        > GDC_bowtie1_index.alt/GRCh38.d1.vd1.fa

Then we build the Bowtie1 index by:

    cd /diskmnt/Projects/cptac_scratch
    mkdir GDC_bowtie1_index.alt
    bowtie-build --threads 8 --seed 201903 \
        GDC_bowtie1_index.alt/GRCh38.d1.vd1.fa \
        GDC_bowtie1_index.alt/GRCh38_d1_vd1_bowtie1_index \
        2> GDC_bowtie1_index.alt/build_index.log 1>&2

Run MapSplice alignment on all samples:

    snakemake mapsplice_all_samples
