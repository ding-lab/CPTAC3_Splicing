#!/bin/bash
set -euo pipefail

readonly GENOME_FA='/diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa'
readonly PER_CHROM_FOLDER='/diskmnt/Projects/cptac_scratch/GRCh38.d1.vd1_per_chrom'

mkdir -p $PER_CHROM_FOLDER
cd $PER_CHROM_FOLDER
csplit -s -z $GENOME_FA '/>/' '{*}'
for csplit_filename in xx*; do
    chrom_name=$(head -n1 $csplit_filename | sed 's/^>//; s/ .*$//')
    mv $csplit_filename "$chrom_name.fa"
done
