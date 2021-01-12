#!/usr/bin/env python3 
# -*- coding:utf-8 -*-
# @Time  : 2020/11/23 13:40
# @Author: lify
# @File  : SAparse.py


import pysam as ps
import networkx as nx
import time
from itertools import groupby
import datetime
from scipy import stats
import numpy as np
from collections import defaultdict
from progressbar import *
from queue import Queue



class Segment(object):
    '''
    Modified from https://github.com/brentp/cigar
    '''
    def __init__(self, pos, cigar):
        self.ref_start = int(pos) - 1
        self.ref_end = self.ref_start
        read_consuming_ops = ("M", "I")
        ref_consuming_ops = ("M", "D")
        cig_iter = groupby(cigar, lambda c: c.isdigit())
        self.read_start, self.read_end = 0, 0
        for i, (g, n) in enumerate(cig_iter):
            counts, tag = int("".join(n)), "".join(next(cig_iter)[1])
            if i == 0 and tag == 'S':
                self.read_start += counts
                self.read_end += counts
            if tag in read_consuming_ops:
                self.read_end += counts
            if tag in ref_consuming_ops:
                self.ref_end += counts

def longest_path(G):
    dist = {} # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0]+1,v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node,(length,_)  = max(dist.items(), key=lambda x:x[1])
    path = []
    while length > 0:
        path.append(node)
        length,node = dist[node]
    return list(reversed(path))

def get_continuous_SAinfo(ont_bam,read_gap):
    '''
    Parse SAinfo from minimap2 aligner of circle-seq nanopore data to find continuous SA for junction detection
    '''
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Parse SAinfo for continuous SA reads\n")
    bamFile = ps.AlignmentFile("%s" % ont_bam, "rb")
    ForwardSeq = {}
    SAInfo = {}
    LinearReads = []
    for read in bamFile:
        if read.is_unmapped:  # unmapped reads
            continue
        if read.is_supplementary:  # supplementary reads
            continue
        readid = read.query_name
        chr1 = bamFile.get_reference_name(read.reference_id)
        strand1 = '+' if not read.is_reverse else '-'
        if not read.has_tag('SA'):
            LinearReads.append(readid)
            if read.get_forward_sequence == None :
                continue
            else :
                ForwardSeq[readid] = read.get_forward_sequence()
        else :
            saInfo = read.get_tag('SA').split(';')[:-1]
            loc = [read.query_alignment_start, read.query_alignment_end, chr1, strand1,
                   read.reference_start, read.reference_end]
            segments = [loc]
            for sa in saInfo:
                chr2, pos, strand2, cigar = sa.split(',')[:4]
                segment = Segment(pos=pos, cigar=cigar)
                segments.append([segment.read_start, segment.read_end, chr2, strand2,
                                 segment.ref_start, segment.ref_end])
            segments.sort()
            SAgraph = nx.DiGraph()
            SAgraph.add_node('END')
            for item in segments :
                nodeID = '-'.join([str(item[0]),str(item[1])])
                SAgraph.add_node(nodeID,segS = item[0],segE = item[1],chrom = item[2], strand = item[3],refS = item[4],refE = item[5],length = item[1]-item[0])
                SAgraph.add_edge(nodeID,'END',weight = (item[1]-item[0]))
            checkedseg = []
            for item in segments:
                for item2 in checkedseg :
                    if ( abs(item2[1] - item[0]) <= read_gap) and item[1] > item2[1] :
                        node1 = '-'.join([str(item[0]), str(item[1])])
                        node2 = '-'.join([str(item2[0]), str(item2[1])])
                        SAgraph.add_edge(node2, node1,weight = (item2[1]-item2[0]))
                checkedseg.append(item)
            LPath = nx.dag_longest_path(SAgraph,weight='weight')
            if len(LPath) == 2 :
                LinearReads.append(readid)
                ForwardSeq[readid] = read.get_forward_sequence()
            elif len(LPath) > 2 :
                info = []
                alignedlen = 0
                for item in LPath :
                    if item != 'END' :
                        alignedlen += SAgraph.nodes[item]['length']
                        info.append([SAgraph.nodes[item]['segS'], SAgraph.nodes[item]['segE'],SAgraph.nodes[item]['chrom'],SAgraph.nodes[item]['strand'],
                                    SAgraph.nodes[item]['refS'],SAgraph.nodes[item]['refE'],SAgraph.nodes[item]['length']])
                forwardSeq = read.get_forward_sequence()
                readlen = read.query_length
                ratio = '%.2f' % (alignedlen/readlen)
                ForwardSeq[readid] = forwardSeq
                SAInfo[readid] = info
#                outfile.write('\t'.join([readid, ratio, str(len(info)), str(info), forwardSeq]) + '\n')
#   ForwardSeq is useless
    return(SAInfo, LinearReads)






class AlignSeg(object):
    '''
    aligned segment information : [readS, readE, chrom, strand, refS, refE, length]
    '''
    def __init__(self, info):
        self.readS = info[0]
        self.readE = info[1]
        self.chrom = info[2]
        self.strand = info[3]
        self.refS = info[4]
        self.refE = info[5]
        self.len = info[6]

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def dfs(DG, nodei, color, circlenode, revisitnode) :
    color[nodei] = 'grey'
    circlenode.append(nodei)
#    print('visiting ' + nodei)
    for nodej in DG.adj[nodei] :
        if color[nodej] == 'grey':
#            print('grey back to ' + nodej)
            revisitnode.append(nodej)
            return True
        if color[nodej] == 'white' and dfs(DG, nodej, color, circlenode, revisitnode) == True:
            return True
    color[nodei] = 'black'
    circlenode.pop()
#    print(nodei + ' is black')
    return False

def hascircle(DG,circlenode, revisitnode) :
    """
        color = white hasn't been visited
        color = grey be visiting
        color = black has been dfs, all children of the spot has been visited
        DG : networkx object
        circlenode = [] all visiting nodes
        revisitnode = [] the node which was revisited
    """
    color = {}
    for i in DG.adj :
        color[i] = 'white'
#    print(color)
    for i in DG.adj :
        if color[i] == 'white' :
            if dfs(DG, i, color, circlenode, revisitnode) == True :
                return True
    return False
#circlenode = []
#revisitnode = []
#hascircle(DG,circlenode, revisitnode)

def parse_junction_from_SAinfo(SAInfo, merge_dist) :
    '''
    Parse junction from continuous SA for junction detection to Tag_Full & Tag_ToAssembly(BreakBSJ + FSJ + interchrom + Unclassified)
    '''
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Parse continuous SAinfo for junction detection\n")
    Tag = {}  # FullBSJ, BreakBSJ, FSJ, Interchrom, Unclassified, Diffstrand, Mixed
    '''
    FullBSJ        # Reads with Full sequence between BSJ
    BreakBSJ       # 2 SA Reads with BSJ but without full sequence
    FSJ            # 2 SA Reads with FSJ
    Interchrom     # 2 SA reads with interchrom junction -> incorrect junction by alignment
    Unclassified   # 2 SA reads with unknown situation
    Diffstrand     # 2 SA reads with same chrom but different strand -> incorrect junction by alignment
    Mixed          # >=3 SA reads without BSJ
    FullJunc       # >=3 SA reads with Full junction
    '''
    Fullinfo = {}  # Junction info of Full type
    BreakBSJ = {}
    for readid in SAInfo.keys() :
        if len(SAInfo[readid]) < 2 :
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),"An error happenend during execution. Exiting")
            sys.exit()
        elif len(SAInfo[readid]) == 2 :
            segment1 = AlignSeg(SAInfo[readid][0])
            segment2 = AlignSeg(SAInfo[readid][1])
            if segment1.chrom != segment2.chrom :
                Tag[readid] = 'Interchrom'
            else :
                if segment1.strand != segment2.strand :
                    Tag[readid] = 'Diffstrand'
                else :
                    if segment1.refS > segment2.refE :
                        Tag[readid] = 'BreakBSJ'
                        BreakBSJ[readid] = ','.join([segment2.chrom, str(segment2.refS), str(segment1.refE), segment2.strand])
                    elif segment1.refE > segment2.refS and segment1.refS <= segment2.refE :
                        Fullinfo[readid] = [','.join(['S', segment2.chrom, segment2.strand, str(segment2.refS)]), ','.join(['E', segment1.chrom, segment1.strand, str(segment1.refE)])]
                        Tag[readid] = 'FullBSJ'
                    elif segment1.refE < segment2.refS :
                        Tag[readid] = 'FSJ'
                    else :
                        Tag[readid] = 'Unclassified'
        else :
            Start = {}
            End = {}
            nodeS = {}
            nodeE = {}
            info = SAInfo[readid]
            SAnum = len(info)
            End[','.join([info[0][2],info[0][3]])] = [info[0][5]]
            Start[','.join([info[SAnum-1][2],info[SAnum-1][3]])] = [info[SAnum-1][4]]
            for subscript in range(1, SAnum - 1):
                seginfo = info[subscript]
                if ','.join([seginfo[2],seginfo[3]]) in Start :
                    Start[','.join([seginfo[2],seginfo[3]])].append(seginfo[4])
                else :
                    Start[','.join([seginfo[2], seginfo[3]])] = [seginfo[4]]
                if ','.join([seginfo[2],seginfo[3]]) in End :
                    End[','.join([seginfo[2],seginfo[3]])].append(seginfo[5])
                else :
                    End[','.join([seginfo[2], seginfo[3]])] = [seginfo[5]]
            for chrom in Start :
                clusters = cluster(Start[chrom],merge_dist)
                for clusteri in clusters :
                    modei = stats.mode(clusteri)[0][0]
                    for starti in clusteri :
                        nodeS[','.join([chrom,str(starti)])] = ','.join(['S',chrom,str(modei)])
            for chrom in End :
                clusters = cluster(End[chrom],merge_dist)
                for clusteri in clusters :
                    modei = stats.mode(clusteri)[0][0]
                    for starti in clusteri :
                        nodeE[','.join([chrom,str(starti)])] = ','.join(['E',chrom,str(modei)])
            DG = nx.DiGraph()
            lastnode = nodeE[','.join([info[0][2],info[0][3],str(info[0][5])])]
            for subscript in range(1, SAnum - 1):
                seginfo = info[subscript]
                nodestart = nodeS[','.join([seginfo[2],seginfo[3],str(seginfo[4])])]
                nodeend = nodeE[','.join([seginfo[2],seginfo[3],str(seginfo[5])])]
                DG.add_edge(lastnode,nodestart,weight=1)
                DG.add_edge(nodestart,nodeend,weight=1000)
                lastnode = nodeend
            visitingnode = []
            revisitnode = []
            if hascircle(DG, visitingnode, revisitnode) == True :
                circlenode = False
                Tag[readid] = 'FullJunc'
                for node in visitingnode :
                    if [node] == revisitnode :
                        circlenode = True
                        Fullinfo[readid] = [node]
                    elif circlenode == True :
                        Fullinfo[readid].append(node)
                    else :
                        continue
            else :
                Tag[readid] = 'Mixed'
    return(Fullinfo, Tag, BreakBSJ)






def segments_cluster(data, maxgap):
    '''
    Arrange data into groups where successive elements
    differ by no more than *maxgap*
    '''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        newGroup = True
        for group in groups[-5:] :
            if (x[0] == group[-1][0]) and (abs(x[1] - group[-1][1]) <= maxgap) and (abs(x[2] - group[-1][2]) <= maxgap) :
                group.append(x)
                newGroup = False
        if newGroup == True :
            groups.append([x])
    return groups

def FullJuncMerge(Fullinfo, Tag, merge_dist) :
    '''
    :param Fullinfo:
    :param SAInfo:
    :param Tag:
    :return:
    No need to take strand information into account while merging
    '''
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Merge and Count Tag Full junctions\n")
    onesegments = []
    multiseg = defaultdict(list)
    for readid in Fullinfo :
        if len(Fullinfo[readid]) == 2 :
            if Fullinfo[readid][0].split(',')[0] == 'S' :
                segment = [Fullinfo[readid][0].split(',')[1],int(Fullinfo[readid][0].split(',')[3]),
                           int(Fullinfo[readid][1].split(',')[3]), readid, Fullinfo[readid][0].split(',')[2]]
            elif Fullinfo[readid][0].split(',')[0] == 'E' :
                segment = [Fullinfo[readid][1].split(',')[1],int(Fullinfo[readid][1].split(',')[3]),
                           int(Fullinfo[readid][0].split(',')[3]), readid, Fullinfo[readid][1].split(',')[2]]
            else :
                Tag[readid] = 'Unclassified'
                continue
            onesegments.append(segment)
        elif len(Fullinfo[readid])%2 == 0 :
            segmentNum = int(len(Fullinfo[readid])/2)
            segments = []
            if Fullinfo[readid][0].split(',')[0] == 'S':
                for st in range(0,segmentNum) :
                    segmenti = []
                    segmenti.append(Fullinfo[readid][2*st].split(',')[1])
                    segmenti.append(Fullinfo[readid][2 * st].split(',')[3])
                    segmenti.append(Fullinfo[readid][2 * st + 1].split(',')[3])
                    segmenti.append(Fullinfo[readid][2 * st + 1].split(',')[2])
                    segments.append(segmenti)
            elif Fullinfo[readid][0].split(',')[0] == 'E':
                if (Fullinfo[readid][0].split(',')[1] == Fullinfo[readid][-1].split(',')[1]) and (Fullinfo[readid][0].split(',')[2] == Fullinfo[readid][-1].split(',')[2]) :
                    segmenti = []
                    segmenti.append(Fullinfo[readid][-1].split(',')[1])
                    segmenti.append(Fullinfo[readid][-1].split(',')[3])
                    segmenti.append(Fullinfo[readid][0].split(',')[3])
                    segmenti.append(Fullinfo[readid][0].split(',')[2])
                    segments.append(segmenti)
                    for st in range(0, segmentNum-1):
                        segmenti = []
                        segmenti.append(Fullinfo[readid][2 * st + 1].split(',')[1])
                        segmenti.append(Fullinfo[readid][2 * st + 1].split(',')[3])
                        segmenti.append(Fullinfo[readid][2 * st + 2].split(',')[3])
                        segmenti.append(Fullinfo[readid][2 * st + 2].split(',')[2])
                        segments.append(segmenti)
            else:
                Tag[readid] = 'Mixed'
                continue
            if segmentNum == 2 :
                segments.sort()
            segments.append([readid])
            multiseg[segmentNum].append(segments)
        else :
            Tag[readid] = 'Mixed'
            continue
    onesegGroup_iteration1 = segments_cluster(onesegments, merge_dist)
    readslist = {}
    for clusteri in onesegGroup_iteration1 :
        chrom = clusteri[0][0]
        starts = []
        ends = []
        reads = []
        strands = []
        for seg in clusteri :
            starts.append(seg[1])
            ends.append(seg[2])
            reads.append(seg[3])
            strands.append(seg[4])
        StartPos = stats.mode(starts)[0][0]
        EndPos = stats.mode(ends)[0][0]
        Strand = stats.mode(strands)[0][0]
        junc = '\t'.join([chrom,  str(StartPos), str(EndPos), Strand])
        readslist[junc] = reads
    return( readslist, multiseg)


def Coverage_count(ont_bam, Junctions) :
    '''
    :param ont_bam:
    :param Junctions:
    :return:
    '''
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Coverage count around junctions\n")
    bamFile = ps.AlignmentFile("%s" % ont_bam, "rb")
    reference_contigs = bamFile.header['SQ']
    header_dict = {}
    for reference in reference_contigs:
        header_dict[reference['SN']] = reference['LN']
    ratiocovL = {}
    ratiocovR = {}
    mannwhitneyuL = {}
    mannwhitneyuR = {}
    cov_range = 50
    num = 0
    total = len(Junctions)
    pbar = ProgressBar().start()
    if 'chrM' in header_dict:
        covchrm = bamFile.count_coverage(contig='chrM')
        sumchrM = list(np.uint32(np.array([covchrm[0], covchrm[1], covchrm[2], covchrm[3]]).sum(axis=0)))
    for junction in Junctions :
        chrom = junction.split('\t')[0]
        junL = int(junction.split('\t')[1])
        junR = int(junction.split('\t')[2])
        if chrom == 'chrM' :
            sumLout = sumchrM[max(junL - cov_range, 0):junL]
            sumLout = [0] * (cov_range - len(sumLout)) + sumLout
            sumLin = sumchrM[junL:min(junL+cov_range, header_dict[chrom])]
            sumLin = sumLin + [0] * (cov_range - len(sumLin))
            sumRout = sumchrM[junR:min(junR+cov_range, header_dict[chrom])]
            sumRout = sumRout + [0] * (cov_range - len(sumRout))
            sumRin = sumchrM[max(junR-cov_range, 0):junR]
            sumRin = [0] * (cov_range - len(sumRin)) + sumRin
        else :
            if junL == 0 :
                sumLout = [0]*cov_range
            else :
                covLout = bamFile.count_coverage(contig=chrom, start=max(junL - cov_range, 0), stop=junL)
                sumLout = list(np.uint32(np.array([covLout[0], covLout[1], covLout[2], covLout[3]]).sum(axis=0)))
                sumLout = [0]*(cov_range-len(sumLout)) + sumLout
            covLin = bamFile.count_coverage(contig=chrom, start=junL, stop=min(junL+cov_range, header_dict[chrom]))
            sumLin = list(np.uint32(np.array([covLin[0], covLin[1], covLin[2], covLin[3]]).sum(axis=0)))
            sumLin = sumLin + [0]*(cov_range-len(sumLin))
            if junR == header_dict[chrom] :
                sumRout = [0]*cov_range
            else :
                covRout = bamFile.count_coverage(contig=chrom, start=junR, stop=min(junR+cov_range, header_dict[chrom]))
                sumRout = list(np.uint32(np.array([covRout[0], covRout[1], covRout[2], covRout[3]]).sum(axis=0)))
                sumRout = sumRout + [0] * (cov_range - len(sumRout))
            covRin = bamFile.count_coverage(contig=chrom, start=max(junR-cov_range, 0), stop=junR)
            sumRin = list(np.uint32(np.array([covRin[0], covRin[1], covRin[2], covRin[3]]).sum(axis=0)))
            sumRin = [0]*(cov_range-len(sumRin)) + sumRin
        ratioL = sum(sumLin) / (sum(sumLout) + sum(sumLin) + 1)
        ratiocovL[junction] = ratioL
        if sumLout == sumLin:
            mannwhitneyuL[junction] = 1
        else:
            mannwhitneyuL[junction] = stats.mannwhitneyu(sumLin, sumLout, alternative='greater').pvalue
        ratioR = sum(sumRin) / (sum(sumRout) + sum(sumRin) + 1)
        ratiocovR[junction] = ratioR
        if sumRout == sumRin:
            mannwhitneyuR[junction] = 1
        else:
            mannwhitneyuR[junction] = stats.mannwhitneyu(sumRin, sumRout, alternative='greater').pvalue
        num += 1
        pbar.update(int((num / (total - 1)) * 100))
        #print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), junction + '\t' + str(num))
    return(ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR)


def Fullinfo_write(Fulloutfile, Junctions, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, Label) :
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Write Single Segments junctions to OnesegFullJunction.out\n")
    FulloutFile = open(Fulloutfile,"w+")
    label = Label + '_Full'
    for junction in Junctions :
        record = '\t'.join([junction, label, str(len(Junctions[junction])), str(mannwhitneyuL[junction]), str(mannwhitneyuR[junction]),
                            str(ratiocovL[junction]), str(ratiocovR[junction]), str(Junctions[junction])])
        FulloutFile.write(record + '\n')
    FulloutFile.close()
    return True





def reads_for_assembly(outfile, LinearReads, Tag) :
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Write reads list for assembly\n")
    outFile = open(outfile,"w+")
    for readid in LinearReads :
        outFile.write(readid + '\n')
    for readid in Tag :
        if Tag[readid] == 'FullJunc' :
            continue
        elif Tag[readid] == 'FullBSJ' :
            continue
        elif Tag[readid] == 'Diffstrand' :
            continue
        elif Tag[readid] == 'Unclassified' :
            continue
        else :
            outFile.write(readid + '\n')
    outFile.close()
    return True





def multisegs_cluster(data, maxgap) :
    data.sort()
    groups = [[data[0]]]
    segnum = len(data[0])
    for x in data[1:]:
        newGroup = True
        for group in groups[-10:] :
            inGroup = True
            for seg in range(segnum) :
                if (x[seg][0] == group[-1][seg][0]) and (abs(x[seg][1] - group[-1][seg][1]) <= maxgap) and \
                        (abs(x[seg][2] - group[-1][seg][2]) <= maxgap) :
                    continue
                else :
                    inGroup = False
                    break
            if inGroup == True :
                newGroup = False
                group.append(x)
        if newGroup == True :
            groups.append([x])
    return groups

def Multiple_Segment_Merge(multiseg, merge_dist) :
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Merge and Count Multiple Segments junctions\n")
    if len(multiseg) == 0 :
        return True
    else :
        Seginfo = {}
        Segreads = defaultdict(list)
        for segnum in multiseg :
            if len(multiseg[segnum]) < 2 :
                continue
            else :
                allsegs = []
                segreads = defaultdict(list)
                for segments in multiseg[segnum] :
                    minseg = min(segments[0:segnum])
                    segqueue = Queue(segnum)
                    for minnum in range(0,segnum) :
                        if segments[minnum] == minseg :
                            break
                        else :
                            segqueue.put(segments[minnum])
                    sortedSegs = segments[minnum:segnum]
                    for item in range(0,segqueue.qsize()) :
                        sortedSegs.append(segqueue.get())
                    for segment in sortedSegs :
                        segment[1] = int(segment[1])
                        segment[2] = int(segment[2])
                    segreads[str(sortedSegs)].append(segments[segnum])
                    allsegs.append(sortedSegs)
                segsgroup = multisegs_cluster(allsegs, merge_dist)
                for group in segsgroup :
                    Seginfo[str(group[0])] = group
                    for item in group :
                        Key = str(group[0])
                        Segreads[Key].append(segreads[str(item)])
        for junction in Segreads :
            reads = str(Segreads[junction]).replace('[','').replace(']','').replace("'",'').split(', ')
            Segreads[junction] = list(np.unique(reads))
    return(Seginfo, Segreads)


def Multiple_Segment_Write(MulsegFullJunction, Seginfo, Segreads, Label):
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Write Multiple Segments junctions to MulsegFullJunction.out\n")
    FulloutFile = open(MulsegFullJunction, "w+")
    label = Label + '_MultiSeg'
    for junction in Seginfo :
        record = '\t'.join([label, str(len(junction.split('], ['))), junction, str(len(Segreads[junction])), str(Segreads[junction])])
        FulloutFile.write(record + '\n')
    FulloutFile.close()
    return True





def BreakJuncMerge(BreakBSJ, merge_dist):
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Merge and Count Tag Break junctions\n")
    segments = []
    for readid in BreakBSJ:
        info = BreakBSJ[readid].split(',')
        segment = [info[0], int(info[1]), int(info[2]), info[3], readid]
        segments.append(segment)
    segGroup = segments_cluster(segments, merge_dist)
    BreakJunc = {}
    for clusteri in segGroup:
        chrom = clusteri[0][0]
        starts = []
        ends = []
        reads = []
        strands = []
        for seg in clusteri:
            starts.append(seg[1])
            ends.append(seg[2])
            reads.append(seg[4])
            strands.append(seg[3])
        StartPos = stats.mode(starts)[0][0]
        EndPos = stats.mode(ends)[0][0]
        Strand = stats.mode(strands)[0][0]
        junc = '\t'.join([chrom,  str(StartPos), str(EndPos), Strand])
        BreakJunc[junc] = reads
    return BreakJunc


def BreakJunc_Write(OnesegBreakJunction, BreakJunc, Label):
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Write Break junctions to OnesegBreakJunction.out\n")
    BreakoutFile = open(OnesegBreakJunction, "w+")
    label = Label + '_Break'
    for junction in BreakJunc :
        record = '\t'.join([junction, label, str(len(BreakJunc[junction])), str(BreakJunc[junction])])
        BreakoutFile.write(record + '\n')
    BreakoutFile.close()
    return True





























