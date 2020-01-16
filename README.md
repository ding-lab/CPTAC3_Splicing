# CPTAC3 RNA-seq splicing/transcript pipeline

version 1.2
	Y2.b3 was run on version 1.1; currently updated to version 1.2
	in version 1.2,  bam files are automatically remove once the analysis for a given sample is complete; output is not affected

This pipeline finds all the transcripts (splice junctions) expressed in each sample. It runs the following set of tools:

- STAR + StringTie



## Setup
Create all the necessary packages inside a conda environment:

    conda create -n cptac3_splice \
        python=3.7 \
        snakemake-minimal \
        star=2.6.1d stringtie=1.3 \
        samtools=1.9 htslib=1.9 \
        gffcompare=0.10



Prepare a list of cases to run the pipeline. The pipeline will run all the
samples belong to the listed case. Then set the `CASE_LIST_PTH` variable inside
`Snakefile`.  For example, `case.list` in year2 batch4 is created by:

    # 	cat /home/mwyczalk_test/Projects/CPTAC3/CPTAC3.catalog/CPTAC3.cases.dat |grep "Y2.b4-RNA"|cut -f 1 > case.list



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

StringTie will produce the following transcript GTFs:
- Per sample at `processed_data/stringtie/{sample}/transcripts.gtf`
- All samples merged at `processed_data/stringtie/merged.gtf`



## Annotations
The pipeline uses GDC hg38 genome reference `GRCh38.d1.vd1`.

STAR, StringTie, and MapSplice all use the same transcript annotation, GENCODE v29 comprehensive annotations on reference chromosomes only (CHR) ([GTF link][gencode-gtf]).

[gencode-gtf]: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz



## Version Change Catalogs
version 1.0 - 1.2 all utlizes star and stringtie for the final output; version changes for the pipeline does not affect the the result of the final output

In version 1.0, processed data are organized by tools (star, stringtie), such that the folder structure is ./processed_data/${tool_name}/${sample_name}_${tissue_type}
In version 1.1, processed data are organized by samples, such that the folder structure is ./processed_data/${sample_name}_${tissue_type}/${tool_name}
In version 1.2, processed data are organized by samples, such that the folder struccutre is ./processed_data/${sample_name}_${tissue_type}/${tool_name}; all the intermediate bam files are removed after job is done so as to free up the space:wq

