#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 04 Dec 2017 00:07:59

import os,sys,argparse

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--refgene',type=str,dest="refgene",help="refgene table downloaded from UCSC genome browser")
    p.add_argument('--hk_gene',type=str,dest="hk_gene",help="housekeeping genes downloaded from http://www.tau.ac.il/~elieis/HKG")
    p.add_argument('--output',type=str,dest="output",help="output file name prefix")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

class Gene(object):
    def __init__(self, gene_name, gene_id):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.chrom = None
        self.start = None
        self.stop = None
        self.strand = None
    def update_pos(self, chrom, start, stop, strand):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand

def parse_hk_gene(filename):
    hk_gene_table = {}
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            row = line.strip().split()
            gene_name = row[0]
            gene_id = row[1]
            hk_gene_table[gene_id] = Gene(gene_name, gene_id)
    fin.close()
    return hk_gene_table

def parse_refGene(filename):
    refgene_list = []
    num = 0
    header = {}
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            row = line.strip().split('\t')
            if num == 0:
                for nn in range(len(row)):
                    header[row[nn]] = nn
            else:
                gene_name = row[header['name2']]
                gene_id = row[header['name']]
                chrom = row[header['chrom']]
                start = int(row[header['txStart']])
                stop = int(row[header['txEnd']])
                strand = row[header['strand']]
                cds_status = row[header['cdsStartStat']]
                if cds_status != 'cmpl':
                    continue
                refgene = Gene(gene_name, gene_id)
                refgene.update_pos(chrom,start,stop,strand)
                refgene_list.append(refgene)
            num += 1
    fin.close()
    return refgene_list

def get_non_hk_gene(hk_gene_table, refgene_list):
    '''
    select genes from refgene_list that are not in hk_gene_list
    1. we will skip genes if either its gene_name or gene_id matches with gene_name or gene_id from hk_gene_list
    2. update hk_genes' genomic position (if gene_id cannot be found in refgene_list, that hk gene will be skipped)
    3. in the case multiple transcripts exist for one gene_name, we will only keep the longest transcript for that gene
    '''
    non_hk_gene_table = {}
    # a dictionary with key as housekeeping gene_name
    hk_gene_name = dict((hk_gene_table[gene_id].gene_name, True) for gene_id in hk_gene_table.keys())
    TOTAL = 0
    MATCH_id = 0
    MATCH_name = 0
    for gene in refgene_list:
        TOTAL += 1
        # if gene_id in hk_gene
        try:
            hk_gene = hk_gene_table[gene.gene_id]
            # update gene position
            hk_gene_table[gene.gene_id].update_pos(gene.chrom, gene.start, gene.stop, gene.strand)
            MATCH_id += 1
            continue
        except KeyError:
            pass
        # if gene_name in hk_gene
        try:
            assert hk_gene_name[gene.gene_name] 
            MATCH_name += 1
            continue
        except KeyError:
            pass
        # the rest is non_hk_gene
        try:
            gene_size = gene.stop - gene.start
            if gene_size > (non_hk_gene_table[gene.gene_name].stop - non_hk_gene_table[gene.gene_name].start):
                non_hk_gene_table[gene.gene_name] = gene
            else:
                pass
        except KeyError:
            non_hk_gene_table[gene.gene_name] = gene
    # update hk_gene_list
    hk_gene_list = []
    for gene_id in hk_gene_table.keys():
        gene = hk_gene_table[gene_id]
        if gene.chrom is not None:
            hk_gene_list.append(gene)
    # get non-hk_gene_list
    non_hk_gene_list = sorted(non_hk_gene_table.keys())
    # log
    print >>sys.stderr, "%d housekeeping genes are selected from %d raw housekeeping genes." % (len(hk_gene_list), len(hk_gene_table))
    print >>sys.stderr, "%d non-housekeeping genes are selected from %d refSeq gene. %d refgene are skipped because of same gene_id as hk genes. %d refgene are skipped because of same gene_name as hk genes." % (len(non_hk_gene_list), len(refgene_list), MATCH_id, MATCH_name)
    return hk_gene_list, non_hk_gene_list

def main():
    global args
    args = parse_arg()
    # parse housekeeping genes
    hk_gene_table = parse_hk_gene(args.hk_gene)
    # parse refGene
    refgene_list = parse_refGene(args.refgene)
    # filter
    hk_gene_list, non_hk_gene_list = get_non_hk_gene(hk_gene_table, refgene_list)
    # report
    with open(args.output+'_hk.bed', 'w') as fout:
        for gene in hk_gene_list:
            print >>fout, gene
    fout.close()
    with open(args.output+'_non_hk.bed', 'w') as fout:
        for gene in non_hk_gene_list:
            print >>fout, gene
    fout.close()
 
if __name__=="__main__":
    main()
