#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# @Time  : 2020/12/15 21:02
# @Author: lify
# @File  : junction_comparing.py


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

def junction_comparing(Abed, Bbed, merge_dist):
    chrom_junctions = defaultdict(list)
    ATag = {}
    BTag = {}
    with open(Abed, 'r') as JuncAFile :
        for line in JuncAFile :
            A = Junction(line.strip().split('\t'))
            ATag[A] = 'Uniq'
            chrom_junctions[A.chrom].append(A)
    with open(Bbed, 'r') as JuncBFile:
        for line in JuncBFile:
            B = Junction(line.strip().split('\t'))
            BTag[B] = 'Uniq'
            for A in chrom_junctions[B.chrom] :
                if (abs(A.start - B.start) <= merge_dist) and (abs(A.end - B.end) <= merge_dist) :
                    BTag[B] = 'Same'
                    ATag[A] = 'Same'
    return(ATag, BTag)

def write_into_file(ATag, BTag, Alabel, Blabel) :
    sameA = open(''.join([Alabel, '_', Blabel, '.same']),"w+")
    sameB = open(''.join([Blabel, '_', Alabel, '.same']),"w+")
    uniqA = open(''.join([Alabel, '_', Blabel, '.uniq']),"w+")
    uniqB = open(''.join([Blabel, '_', Alabel, '.uniq']), "w+")
    for item in ATag :
        junction = '\t'.join([item.chrom, str(item.start), str(item.end)])
        if ATag[item] == 'Same' :
            sameA.write(junction + '\n')
        elif ATag[item] == 'Uniq' :
            uniqA.write(junction + '\n')
        else :
            continue
    for item in BTag :
        junction = '\t'.join([item.chrom, str(item.start), str(item.end)])
        if BTag[item] == 'Same' :
            sameB.write(junction + '\n')
        elif BTag[item] == 'Uniq' :
            uniqB.write(junction + '\n')
        else :
            continue
    sameA.close()
    sameB.close()
    uniqA.close()
    uniqB.close()
    return True



ATag, BTag = junction_comparing(FileA, FileB, merge_dist)
write_into_file(ATag, BTag, labelA, labelB)