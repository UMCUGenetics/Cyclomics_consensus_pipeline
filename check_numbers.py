#! /usr/bin/env python
import sys, os, re
import glob
import subprocess
import gzip

from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
    parser = OptionParser();
    group = OptionGroup(parser, "Main options")
    group.add_option("-r", dest="rawdir", metavar="[PATH]", help="full path to raw data folder (FASTQ location)")
    group.add_option("-p", dest="prodir", metavar="[PATH]", help="full path to processed data folder")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    if opt.rawdir:
        rawdir = opt.rawdir
    else:
        sys.exit("provide folder to raw data" )

    if opt.prodir:
        prodir = opt.prodir
    else:
        sys.exit("provide folder to processed data" )

    files = subprocess.getoutput("find {rawdir}  -iname \"*fastq\"".format(rawdir=rawdir).split()
    if len(files) == 0: ## likely that fastq are compressed
        fastq = subprocess.getoutput("find {rawdir} -iname \"*fastq.gz\" -exec gunzip -cd {} \; | wc -l".format(rawdir=rawdir)
    else:
        fastq = subprocess.getoutput("find {rawdir} -iname \"*fastq\" -exec cat {} \; | wc -l").format(rawdir=rawdir)

    reads = int(int(fastq)/4)
    bams = subprocess.getoutput("find {prodir}/bam/ -iname \"*tar\" -exec tar -tf {} \; |grep bai | wc -l").format(prodir=prodir)
    m5 = subprocess.getoutput("find  {prodir}/m5/ -iname \"*tar\" -exec tar -tf {} \; | wc -l").format(prodir=prodir)
    consensus = subprocess.getoutput("find {prodir}/consensus/ -iname \"*tar\" -exec tar -tf {} \; | wc -l").format(prodir=prodir)


    print("reads in fastq file ", reads)
    print("bai files in bam tars  ", bams)
    print("m5 files in m5 tars  ", m5)
    print("consensus files in consensus tars  ", consensus)
    if int(reads) == int(bams) and  int(reads) == int(m5) and int(reads) == int(consensus) and int(bams) == int(m5) and int(bams) == int(consensus) and  int(m5) == int(consensus):
        print("####\nEverything OK\n####")
    else:
       print("####\nProblem!\n####")

