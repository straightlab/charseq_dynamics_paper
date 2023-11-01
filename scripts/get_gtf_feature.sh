#!/bin/bash

FEATURE="$1"

TRANSCRIPTOME_ROOT="/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/transcriptome_hg38/gencode_annotations"
CHROMS="/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/hg38/hg38.chrom.sizes"

bedtools merge -i <(bedtools slop -i <(cat "$TRANSCRIPTOME_ROOT"/gencode.v28.basic.annotation.gtf | grep $FEATURE) -g "$CHROMS" -b $2)
