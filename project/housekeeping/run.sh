#!/bin/bash

python src/match_gene.py --refgene data/hg19_refGene.txt --hk_gene data/HK_genes.txt --output results/HK
# 3791 housekeeping genes are selected from 3804 raw housekeeping genes.
# 15751 non-housekeeping genes are selected from 52925 refSeq gene. 4028 refgene are skipped because of same gene_id as hk genes. 7561 refgene are skipped because of same gene_name as hk genes.

bedtools makewindows -g data/hg19F.genome -w 20000 >data/hg19F_20k.bed

python Norma/src/TSA-seq_anno.py --bed data/hg19F_20k.bed --anno cell\ type\ compare/HCT116.bw cell\ type\ compare/HFFc6.bw cell\ type\ compare/K562.bw results/HK_hk.bed.gz results/HK_non_hk.bed.gz --mode signal_mean signal_mean signal_mean anno anno --label HCT116 HFF K562 hk non_hk --genome data/hg19.fa --genome_size data/hg19.genome --output results/20kb_anno.txt

python Norma/src/TSA-seq_anno.py --update --bed results/20kb_anno.txt --anno data/hg19_Gap.bed.gz --mode count --label Gap --genome data/hg19.fa --genome_size data/hg19.genome --output results/20kb_anno.txt

Rscript src/housekeeping_gene.R
