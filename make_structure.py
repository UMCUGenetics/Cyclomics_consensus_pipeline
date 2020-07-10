#! /usr/bin/env python
import sys
import os
import re
from optparse import OptionParser, OptionGroup
import glob
import subprocess
import pysam
import numpy as np
import settings

class Segment:
    def __init__( self, id, qname, flag, rname, pos, mapq, length ):
        self.id = id
        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.length = int(length)
        self.end = ( int(pos) + int(length) )
        self.clip = False
        self.pid = False

    def parseCigar( self, cigar ):
        if self.flag & 16:
            self.clip = cigar['H'][1]
        else:
            self.clip = cigar['H'][0]
        self.end += cigar['D']
        self.end -= cigar['I']
        self.pid = format( cigar['=']/float(self.length), '.3f' )

def getKey(item):
    return item[3]

def Phred(basescore):
    return ord(basescore) - 33

def Overlap(list1, list2):
    size_overlap = False
    if list1[0] >= list2[0] and list1[0] <= list2[1]: ## list1 begin falls within list2
        size_overlap = True
    if list1[1] >= list2[0] and list1[1] <= list2[1]: ## list1 stop falls within list2
        size_overlap = True
    return size_overlap

def mad(data, axis=None):
    return np.mean(np.absolute(data - np.median(data, axis)), axis)  ## mean of the deviation from the median

if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("-i", default = "./", dest = "wkdir", help = "full path to folder with BAM files [default = ./")
    group.add_option("-o", default = settings.STRUCTURE_OVERLAP, dest = "overlap", help = "overlap in bp to determine if insert fragements are from the same molecule/amplicon [default STRUCTURE_OVERLAP in settings.py]")
    group.add_option("-f", default = settings.STRUCTURE_FLANK, dest = "flank", help = "flank to determine unmapped regions to either BB or I [default STRUCTURE_FLANK in settings.py")
    group.add_option("-m", default = settings.STRUCTURE_MAD, dest = "mad_threshold", help = " if mad score for insert startsite is more than threshold, report threshold [default STRUCTURE_MAD in settings.py]")
    group.add_option("-r", default = settings.STRUCTURE_INTARGET, dest = "interval", help = "interval region of interest chr:start-stop [default STRUCTURE_INTARGET in settings.py] ")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    overlap = int(opt.overlap)
    flank = int(opt.flank)
    region = re.split("[:-]", opt.interval)
    wkdir = opt.wkdir
    mad_threshold = float(opt.mad_threshold)

    print("ReadName\treadStructure\tMADinsert\tInsertCount\tDifInserts\tBackboneCount\tDifBackbone\tB\tI\tBI\tMeanBaseQualityMappedRead")
    for subfolder in subprocess.getoutput("find " + str(wkdir) + " -mindepth 1 -type d -iname \"*\"").split():
        tarfile = str(subfolder) + "/" + str(subfolder).split("/")[-1] + ".tar"
        if os.path.isfile(tarfile):
            files = subprocess.getoutput("tar -tf " + str(tarfile)).split()
            for f in files:
                if "bai" not in f:
                    os.system("tar -xf " + str(tarfile) + " " + str(f))
                    bamfile = pysam.AlignmentFile(str(f), "rb")
                    rlist = []
                    total = 0
                    length = 0
                    quality = [0]
                    startlist = []
                    for line in bamfile: ## parse BAM file
                        """Determine hardclip start and end, mapped length, original length"""
                        length = line.query_length  	## length in FASTQ
                        length_mapped = line.reference_length ## mapped length
                        if line.flag & 16:
                            if line.cigartuples[-1][0] == 5:
                                start = line.cigartuples[-1][1]
                            else:
                                start = 0
                            if line.cigartuples[0][0] == 5:
                                stop = line.cigartuples[0][1]
                            else:
                                stop = 0
                        else:
                            if line.cigartuples[0][0] == 5:
                                start = line.cigartuples[0][1]
                            else:
                                start = 0
                            if line.cigartuples[-1][0] == 5:
                                stop = line.cigartuples[-1][1]
                            else:
                                stop = 0
                        total = length + stop + start ## total lenth of the sequenced read
                        rlist += [[bamfile.get_reference_name(line.reference_id), line.reference_start, line.flag, start, stop, length, total, length_mapped]]
                        quality += line.query_qualities

                    meanq = np.mean(quality)
                    rlist = sorted(rlist, key=getKey)
                    start = 0
                    stop = total
                    x = 0
                    dic = {"gap":[], "switch":[], "BI":[]}
                    while x < len(rlist):
                        if rlist[x][3] >= start: ## calculate gap between fragments
                            dic["gap"] += [rlist[x][3] - start]
                            start = rlist[x][3] + rlist[x][5] ## make new start end of fragment
                        x += 1
                    dic["gap"] += [stop - start]  #last unmapped region
                    full_sequence = ""
                    location = {}
                    backbone = {}
                    startsite = []
                    z = 0
                    y = 0
                    b = 0
                    i_fragment = [0]
                    b_fragment = [0]
                    def ReturnChunk(length, backbone, ori, chrom, position, mappedl, number):
                         return '{}:{}:{}:{}:{}:{}:{},'.format(length, backbone, ori, chrom, position, mappedl, number)

                    for item in rlist:
                        if dic["gap"][z] > 0: ## print gaps in structure
                            full_sequence += str(dic["gap"][z]) + ":U,"
                        if "BB" in item[0] or "pJet" in item[0] or "pUC" in item[0]:
                            try:
                                backbone[item[0]]
                                b_fragment[backbone[item[0]]] += 1
                            except:
                                backbone[item[0]] = b
                                if b > 0:
                                    b_fragment += [1]
                                else:
                                    b_fragment[backbone[item[0]]] += 1
                                b += 1
                            if item[2] & 16:
                                full_sequence += ReturnChunk(item[5], "BB", "R", item[0], item[1], item[7], backbone[item[0]])
                            else:
                                full_sequence += ReturnChunk(item[5], "BB", "F", item[0], item[1], item[7], backbone[item[0]])
                        else:
                            if str(item[0]) == str(region[0]) and int(item[1])> int(region[1]) and int(item[1]) < int(region[2]): ## if startsite of fragment is within target interval
                                startsite += [int(item[1])]
                            if bool(location) == False:
                                for i in range(int(item[1]) - overlap, int(item[1]) + overlap):
                                    try:
                                         location[i]
                                    except:
                                        location[i] = y
                            try: ## within range, thus considered same amplicon
                                location[item[1]]
                                i_fragment[location[item[1]]] += 1
                                if item[2] & 16:
                                    full_sequence += ReturnChunk(item[5], "I", "R", item[0], item[1], item[7], location[item[1]])
                                else:
                                    full_sequence += ReturnChunk(item[5], "I", "F", item[0], item[1], item[7], location[item[1]])
                            except: ## outside range, thus considered different amplicon
                                y += 1
                                i_fragment += [1]
                                for i in range(int(item[1]) - overlap, int(item[1]) + overlap):
                                    try:
                                        location[i]
                                    except:
                                        location[i] = y
                                if item[2] & 16:
                                    full_sequence += ReturnChunk(item[5], "I", "R", item[0], item[1], item[7], location[item[1]])
                                else:
                                    full_sequence += ReturnChunk(item[5], "I", "F", item[0], item[1], item[7], location[item[1]])
                        z += 1
                    if dic["gap"][z] > 0: ## needed for last unknown sequence
                         full_sequence += str(dic["gap"][z]) + ":U"
                    if dic["gap"][z] == 0: 	## remove last comma (no Unknown insert)
                         full_sequence = full_sequence[0:-1]

                    bb_count = full_sequence.count("BB:") ## ":" is added to prevent additional counting in chromosome name
                    i_count = full_sequence.count("I:")

                    printline = str(f) + "\t"
                    if len(full_sequence)> 0:
                        printline += full_sequence + "\t"
                    else:
                        printline += "nan\t"
                    if len(startsite) > 0:
                        mad_score = "%.2f" % float(mad(startsite))
                        if float(mad_score) <= float(mad_threshold):
                            printline += str("%.2f" % float(mad(startsite))) + "\t"
                        else:
                            printline += str("%.2f" % float(mad_threshold)) + "\t"
                    else:
                        printline += "nan\t"
                    printline += ",".join(map(str, i_fragment)) + "\t"
                    if i_fragment[0] == 0:	## if insert count is zero, report as none
                        printline += "nan" + "\t"
                    else:
                        printline += str(len(i_fragment)) + "\t"

                    printline += ",".join(map(str, b_fragment)) + "\t"
                    if b_fragment[0] == 0:    ## if backbone count is zero, report as none
                        printline += "nan" + "\t"
                    else:
                        printline += str(len(b_fragment)) + "\t"

                    printline += str(bb_count) + "\t"
                    printline += str(i_count) + "\t"
                    try:
                        printline += str("%.3f" % float(float(bb_count)/float(i_count))) + "\t"
                    except:
                        printline += "nan\t"
                    try:
                        printline += str("%.2f" % (meanq))
                    except:
                       printline += "nan"
                    print(printline)
                    os.system("rm " + str(f))
