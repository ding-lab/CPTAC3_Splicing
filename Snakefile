import csv
from collections import namedtuple
from pathlib import Path
import re

MATT_CPTAC3_CATALOG_PTH = '/home/mwyczalk_test/Projects/CPTAC3/CPTAC3.catalog'
SERVER_LOCATION = 'katmai'
CASE_LIST = 'case.list'

GENOME_FA = '/diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa'
STAR_INDEX_FOLDER = '/diskmnt/Datasets/Reference/GDC/star_genome_d1_vd1_gencode_comp_chr_v29'
STAR_GTF = '/diskmnt/Datasets/Reference/GDC/gencode.v29.annotation.gtf'

# Define cases and samples in use {{{
CASES = set(open(CASE_LIST).read().splitlines())
LOCAL_MAP_PTH = Path(MATT_CPTAC3_CATALOG_PTH, f'{SERVER_LOCATION}.BamMap.dat')
# Define Sample object
Sample = namedtuple('Sample', 'case sample_type')

# Select all the available samples of the selected cases
SAMPLES = set()
with open(LOCAL_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        if row['case'] in CASES and row['experimental_strategy'] == 'RNA-Seq':
            s = Sample(row['case'], row['sample_type'])
            SAMPLES.add(s)
# }}}


rule link_cptac_gdc_rna_fastqs:
    """Link the CPTAC RNA-seq FASTQs locally."""
    input: local_map=LOCAL_MAP_PTH
    output: expand('external_data/GDC_RNA_fq/{s.case}_{s.sample_type}.{read}.fastq.gz', \
                   s=SAMPLES, read=['R1', 'R2'])
    run:
        reader = csv.DictReader(open(input['local_map']), dialect='excel-tab')
        for row in reader:
            # Keep only RNA-seq UUIDs
            if not (row['experimental_strategy'] == 'RNA-Seq' and row['data_format'] == 'FASTQ'):
                continue
            s = Sample(row['case'], row['sample_type'])
            if s not in SAMPLES:
                continue
            # Get the RNA-seq read strand
            read = re.search(r'\.(R1|R2)\.\w+$', row['# sample_name']).group(1)
            src_pth = Path(row['data_path'])
            assert src_pth.exists()   # make sure the FASTQ file exists
            dst_pth = Path(f'external_data/GDC_RNA_fq/'
                           f'{s.case}_{s.sample_type}.{read}.fastq.gz')
            dst_pth.symlink_to(Path(src_pth))


# STAR 2-pass alignment {{{
rule star_align_pass1:
    """STAR genome alignemnt pass 1 of one sample."""
    input:
        r1_fq='external_data/GDC_RNA_fq/{sample}.R1.fastq.gz',
        r2_fq='external_data/GDC_RNA_fq/{sample}.R2.fastq.gz'
    output:
        sj='processed_data/star/{sample}/pass1/SJ.out.tab'
    params:
        star_ix=STAR_INDEX_FOLDER,
        star_gtf=STAR_GTF,
        out_folder='processed_data/star/{sample}/pass1/'
    log: 'logs/star_pass1/{sample}.log'
    threads: 4
    shell:
        # Follow the GDC's parameters
        'STAR '
        '--genomeDir {params.star_ix} '
        '--readFilesIn {input.r1_fq} {input.r2_fq} '
        '--readFilesCommand zcat '
        '--runThreadN {threads} '
        '--outFileNamePrefix {params.out_folder} '
        '--outFilterMultimapScoreRange 1 '
        '--outFilterMultimapNmax 20 '
        '--outFilterMismatchNmax 10 '
        '--alignIntronMax 500000 '
        '--alignMatesGapMax 1000000 '
        '--sjdbScore 2 '
        '--alignSJDBoverhangMin 1 '
        '--genomeLoad NoSharedMemory '
        '--outFilterMatchNminOverLread 0.33 '
        '--outFilterScoreMinOverLread 0.33 '
        '--sjdbOverhang 100 '
        '--outSAMstrandField intronMotif '
        '--outSAMtype None '
        '--outSAMmode None '
        '> {log}'


rule star_align_pass2:
    """STAR genome alignemnt pass 2 of one sample."""
    input:
        r1_fq='external_data/GDC_RNA_fq/{sample}.R1.fastq.gz',
        r2_fq='external_data/GDC_RNA_fq/{sample}.R2.fastq.gz',
        all_pass1_sj_tabs=expand(rules.star_align_pass1.output['sj'], \
                                 sample=[f'{s.case}_{s.sample_type}' for s in SAMPLES])
    output:
        bam='processed_data/star/{sample}/pass2/Aligned.sortedByCoord.out.bam'
    params:
        star_ix=STAR_INDEX_FOLDER,
        star_gtf=STAR_GTF,
        out_folder='processed_data/star/{sample}/pass2/'
    log: 'logs/star_pass2/{sample}.log'
    threads: 4
    shell:
        # Run the same command as pass1, but passing the split junction tabs of all samples
        'STAR '
        '--genomeDir {params.star_ix} '
        '--readFilesIn {input.r1_fq} {input.r2_fq} '
        '--readFilesCommand zcat '
        '--runThreadN {threads} '
        '--outFileNamePrefix {params.out_folder} '
        '--outFilterMultimapScoreRange 1 '
        '--outFilterMultimapNmax 20 '
        '--outFilterMismatchNmax 10 '
        '--alignIntronMax 500000 '
        '--alignMatesGapMax 1000000 '
        '--sjdbScore 2 '
        '--alignSJDBoverhangMin 1 '
        '--genomeLoad NoSharedMemory '
        '--outFilterMatchNminOverLread 0.33 '
        '--outFilterScoreMinOverLread 0.33 '
        '--sjdbOverhang 100 '
        # Pass SJ tabs of all samples (main difference to PASS 1)
        '--sjdbFileChrStartEnd {input.all_pass1_sj_tabs} '
        '--limitBAMsortRAM 0 '
        '--outSAMstrandField intronMotif '
        '--outSAMattributes NH HI NM MD AS XS '
        '--outSAMunmapped Within '
        '--outSAMtype BAM SortedByCoordinate '
        # Add SAM version number (likely no effect on anything)
        # Ref: https://www.biostars.org/p/121135/#121138
        '--outSAMheaderHD @HD VN:1.4 '
        '> {log}'

# }}}


rule star_align_all_samples:
    """Run STAR alignment all samples."""
    input: all_bams=expand(rules.star_align_pass2.output['bam'], \
                           sample=[f'{s.case}_{s.sample_type}' for s in SAMPLES])


rule samtools_index_bam:
    """Samtools index one BAM."""
    input: '{name}.bam'
    output: '{name}.bam.bai'
    shell: 'samtools index {input}'
