#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# @Time  : 2020/12/9 16:05
# @Author: lify
# @File  : run.py

import sys
import argparse
sys.path.insert(0, '/mnt/program/fuyu/04_Nanopore/eccDNA_lwx/00_LongReadsCircDNADetection_version5/v5.1')
import SAparse as ont
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("prefix", help="prefix of bam file, \n eg. MCF12A.sorted.bam, then prefix is MCF12A", type=str)
parser.add_argument("Label", help="Label of sample in outfile", type=str)


args = parser.parse_args()
prefix = args.prefix
label = args.Label

ont_bam = prefix + '.sorted.bam'
read_gap = 50
merge_dist = 50
OnesegFullJunction = label + '.DiGraph.OnesegFullJunction.out'
OnesegBreakJunction = label + '.DiGraph.OnesegBreakJunction.out'
MulsegFullJunction = label + '.DiGraph.MulsegFullJunction.out'
reads4assemblefile = label +  '.reads4canu.list'

SAInfo, LinearReads = ont.get_continuous_SAinfo(ont_bam,read_gap)
Fullinfo, Tag, BreakBSJ = ont.parse_junction_from_SAinfo(SAInfo, merge_dist)
Junctions, multiseg = ont.FullJuncMerge(Fullinfo, Tag, merge_dist)

# Full Single Segment Junctions
ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR = ont.Coverage_count(ont_bam, Junctions)
ont.Fullinfo_write(OnesegFullJunction, Junctions, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, label)

# Full Multiple Segments Junctions
Seginfo, Segreads = ont.Multiple_Segment_Merge(multiseg, merge_dist)
ont.Multiple_Segment_Write(MulsegFullJunction, Seginfo, Segreads, label)

# Break Junctions
BreakJunc = ont.BreakJuncMerge(BreakBSJ, merge_dist)
ont.BreakJunc_Write(OnesegBreakJunction, BreakJunc, label)

# Other reads for assembly
ont.reads_for_assembly(reads4assemblefile, LinearReads, Tag)


print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "CircDNA Detection Done!\n")


#seqkit grep --pattern-file MCF12A_2_1_3.reads4canu.list MCF12A_2_1_3.q7.fastq > MCF12A_2_1_3.reads4canu.fastq
#nohup canu -p MCF12A_2_1_3 -d /home/luna/lify/MCF12A_assembly/MCF12A_2_1_3 genomeSize=50M -nanopore-raw MCF12A_2_1_3.reads4canu.fastq &
#awk -F '\t' '{if ($5>1 && $6<=0.01 && $7<=0.01) print $0}' MCF12A_smallData.DiGraph.OnesegFullJunction.out > MCF12A_smallData.DiGraph.OnesegFullJunction.out.2reads.P01













