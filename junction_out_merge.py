#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# @Time  : 2020/12/18 15:32
# @Author: lify
# @File  : junction_out_merge.py

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("P_juncfile", help="primary junction file A", type=str)
parser.add_argument("S_juncfile", help="supplement junction file B", type=str)
parser.add_argument("prefix", help="prefix of outfile", type=str)

args = parser.parse_args()
FileP = args.P_juncfile
FileS = args.S_juncfile
prefix = args.prefix

merge_dist = 100

class Junction(object) :
    def __init__(self, info):
        self.chrom = info[0]
        self.strand = info[3]
        self.start = int(info[1])
        self.end = int(info[2])
        self.label = info[4]
        self.count = int(info[5])
        self.info = '\t'.join(info[6:])



def junction_merging(primary, supplement, merge_dist) :
    chrom_junctions = defaultdict(list)
    Juncount = {}
    SupTag = {}
    with open(primary, 'r') as primaryFile:
        for line in primaryFile:
            priJunc = Junction(line.strip().split('\t'))
            Juncount[priJunc] = priJunc.count
            chrom_junctions[priJunc.chrom].append(priJunc)
    with open(supplement, 'r') as supplementFile:
        for line in supplementFile:
            supJunc = Junction(line.strip().split('\t'))
            SupTag[supJunc] = line.strip()
            for junc in chrom_junctions[supJunc.chrom]:
                if (abs(junc.start - supJunc.start) <= merge_dist) and (abs(junc.end - supJunc.end) <= merge_dist):
                    SupTag[supJunc] = 'Same'
                    Juncount[junc] += 1
    return(Juncount, SupTag)


def write_into_file(Juncount, SupTag, prefix) :
    finalFile = open(''.join([prefix, '.DiGraph.OnesegJunction.merged.out']),"w+")
    for item in Juncount :
        junctioninfo = '\t'.join([item.chrom, str(item.start), str(item.end), item.strand, item.label, str(Juncount[item]), item.info])
        finalFile.write(junctioninfo + '\n')
    for item in SupTag :
        if SupTag[item] != 'Same' :
            junctioninfo = SupTag[item]
            finalFile.write(junctioninfo + '\n')
    finalFile.close()
    return True



Juncount, SupTag = junction_merging(FileP, FileS, merge_dist)
write_into_file(Juncount, SupTag, prefix)


