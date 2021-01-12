#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 2021/1/4 17:33
# @Author: lify
# @File  : junction_intersection.py

# A In B
# A Intersect B
# A Include B
# A Not B
# A Same B

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("A_label", help="label of junction file A", type=str)
parser.add_argument("A_juncfile", help="junction file A", type=str)
parser.add_argument("B_label", help="label of junction file B", type=str)
parser.add_argument("B_juncfile", help="junction file B", type=str)

args = parser.parse_args()
FileA = args.A_juncfile
FileB = args.B_juncfile
labelA = args.A_label
labelB = args.B_label

merge_dist = 100

class Junction(object) :
    def __init__(self, info):
        self.chrom = info[0]
        self.start = int(info[1])
        self.end = int(info[2])

def junction_intersection(Abed, Bbed, merge_dist, labelA, labelB):
    outfile = open(labelA + '_intersect_' + labelB + '.out', "w+")
    chrom_junctions = defaultdict(list)
    with open(Bbed, 'r') as JuncBFile :
        for line in JuncBFile :
            B = Junction(line.strip().split('\t'))
            chrom_junctions[B.chrom].append(B)
    with open(Abed, 'r') as JuncAFile:
        for line in JuncAFile:
            A = Junction(line.strip().split('\t'))
            for B in chrom_junctions[A.chrom] :
                if (abs(A.start - B.start) <= merge_dist) and (abs(A.end - B.end) <= merge_dist) :
                    intersection = '\t'.join([A.chrom, str(A.start), str(A.end), 'Same', B.chrom, str(B.start), str(B.end)])
                    outfile.write(intersection + '\n')
                elif (A.start - B.start >= 0) and (A.end - B.end <= 0):
                    intersection = '\t'.join([A.chrom, str(A.start), str(A.end), 'Included', B.chrom, str(B.start), str(B.end)])
                    outfile.write(intersection + '\n')
                elif (A.start - B.start <= 0) and (A.end - B.end >= 0):
                    intersection = '\t'.join([A.chrom, str(A.start), str(A.end), 'Contain', B.chrom, str(B.start), str(B.end)])
                    outfile.write(intersection + '\n')
                elif (min(A.end, B.end) > max(A.start, B.start)) :
                    intersection = '\t'.join([A.chrom, str(A.start), str(A.end), 'Intersect', B.chrom, str(B.start), str(B.end)])
                    outfile.write(intersection + '\n')
                else:
                    continue
    outfile.close()
    return True

junction_intersection(FileA, FileB, merge_dist, labelA, labelB)