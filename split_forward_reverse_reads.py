#! /usr/bin/env python
import sys
import os
import re
import subprocess
import pysam
from optparse import OptionParser
from optparse import OptionGroup
import settings

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest = "wkdir", metavar = "[PATH]", help = "full path to BAM input folder [default = ./")
        group.add_option("-o", dest = "outdir", metavar = "[PATH]", help = "full path to output folder [default = ./")
        group.add_option("-b", default = settings.bwa, dest = "bwa", metavar = "[PATH]", help = "full path to bwa executable [default bwa in settings.py]")
        group.add_option("-g", dest = "gene", metavar = "[STRING]", help = "gene of interest [default = TP53")
        group.add_option("--mp", default =settings.SPLIT_MIN_PERC, dest = "minperc", metavar = "[FLOAT]", help = "minimum percentage of forward or reverse reads required to include in forward or reverse bin [default SPLIT_MIN_PERC in settings.py")
        group.add_option("--sa", default = settings.sambamba, dest = "sambamba", metavar = "[PATH]", help = "full path to sambamba executable [default sambamba in settings.py]")
        group.add_option("--b5", default = settings.bam2m5, dest = "bam2m5", metavar = "[PATH]", help = "full path to bam2m5 executable [default bam2m5 in settings.py")
        group.add_option("--rf", default = settings.full_ref, dest = "refgenome_full", metavar = "[PATH]", help = "full path to complete reference genome [default full_ref in settings.py]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

infolder = opt.wkdir # should be folder with sorted.bam
outfolder = opt.outdir
minperc = float(opt.minperc)

if outfolder.endswith('/'): # chop last "/" if present
    outfolder = outfolder[0:-1]

if opt.gene:
    gene = str(opt.gene)
else:
    gene = "TP53"

# Log GIT version + commit of repository
if str(sys.argv[0]) == "python":
    repo = "/".join(sys.argv[1].split("/")[0:-1])
else:
    repo = "/".join(sys.argv[0].split("/")[0:-1])

os.system("git --git-dir=" + str(repo) + "/.git describe --tags >" + str(outfolder) + "/GIT_split_forward_reverse_reads.log")
os.system("git --git-dir=" + str(repo) + "/.git log >>" + str(outfolder) + "/GIT_split_forward_reverse_reads.log")

l = ["forward","reverse"]
trim = 0
bwa = str(opt.bwa) + " mem -t 2 -c 100 -M -R"
refgenome_full = opt.refgenome_full
sambamba = opt.sambamba
sample = outfolder.split("/")[-2].split("_")[0]

os.system ("mkdir " + str(outfolder) + "/tmp")
os.system("rm " + str(outfolder) + "/" + str(outfolder.split("/")[-1]) + "_reverse.fasta")
os.system("rm " + str(outfolder) + "/" + str(outfolder.split("/")[-1]) + "_forward.fasta")

for folder in os.listdir(infolder):
    tarfile = str(infolder) + "/" + str(folder) + "/" + str(folder).split("/")[-1] + ".tar"
    if  os.path.isfile(tarfile): # skip empty folders
        files = subprocess.getoutput("tar -tf " + str(tarfile)).split()
        for f in files:
            if "bai" not in f:
                if os.path.isfile(tarfile) == False:
                    break
                os.system("tar -axf " + str(tarfile) + " " + str(f))
                in_file = pysam.AlignmentFile(str(outfolder) + "/" + str(f), "rb")
                forward = 0
                reverse = 0
                total = 0
                for line in in_file:
                    if str(gene) in line.reference_name:
                        if line.flag == 16:
                            reverse += 1
                        else:
                            forward += 1
                        total += 1
                os.system("rm " + str(outfolder) + "/" + str(f))
                conf = f[0:-11] + str(".consensus")
                tarfilec = str(infolder.replace("bam", "consensus")) + "/" + str(folder) + "/" + str(folder).split("/")[-1] + ".tar"
                count_rev_org = 0
                count_for_org = 0
                count_rev_new = 0
                count_for_new = 0
                if reverse > forward:
                    count_rev_org += 1
                elif forward > reverse:
                    count_for_org += 1

                if total > 0: 
                    rev_readperc = (float(reverse)/float(total))*100
                    for_readperc = (float(forward)/float(total))*100
                    if rev_readperc > minperc:
                        os.system("tar -axf " + str(tarfilec) + " " + str(conf) + " -O >> " + str(outfolder.split("/")[-1]) + "_reverse.fasta")
                        count_rev_new += 1
                    elif for_readperc > minperc:
                        os.system("tar -axf " + str(tarfilec) + " " + str(conf) + " -O >> " + str(outfolder.split("/")[-1]) + "_forward.fasta")
                        count_for_new += 1
                    else:
                        pass # discard reads      
                else:
                    rev_readperc = 0
                    for_readperc = 0
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(f, reverse, forward, count_rev_new, count_for_new, rev_readperc, for_readperc))
     
for item in l:
    os.system (str(bwa) + "\"@RG\\tID:"+str(item) + "\\tSM:" + (sample) + "\\tPL:NANOPORE\\tLB:" + str(sample) + "\" " + str(refgenome_full) + " " + str(outfolder.split("/")[-1]) + "_" + str(item) + ".fasta | " + str(sambamba) + " view -S -f bam /dev/stdin | " + str(sambamba) + " sort -t 2 --tmpdir=./tmp /dev/stdin -o " + str(outfolder.split("/")[-1]) + "_" + str(item) + ".sorted.bam")
    os.system("rm " + str(outfolder.split("/")[-1]) + "_" + str(item) + ".fasta")

