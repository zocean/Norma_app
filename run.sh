#!/bin/bash

python src/match_gene.py --refgene data/hg19_refGene.txt --hk_gene data/HK_genes.txt --output results/HK
# 3791 housekeeping genes are selected from 3804 raw housekeeping genes.
# 15751 non-housekeeping genes are selected from 52925 refSeq gene. 4028 refgene are skipped because of same gene_id as hk genes. 7561 refgene are skipped because of same gene_name as hk genes.
