#! /usr/bin/env python
import os
from optparse import OptionParser
from optparse import OptionGroup
import pysam
import settings

if __name__ == "__main__":
    parser = OptionParser();
    group = OptionGroup(parser, "Main options")
    group.add_option("-b", dest = "blacklist", help = "blacklist of reads that should be excluded for basecalling. File should contain 1 column with readID to exclude")
    group.add_option("-s", default = setting.samtools, dest = "samtools", help = "full path to samtools executable [ default samtools in settings.py")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()
   
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith((".sorted.bam")):
            in_file = open(opt.blacklist,"r").readlines()
            output_id = "{bam_id}_bl.bam".format(bam_id = str(f)[0:-4])
            bl_dic = {}
            for line in in_file:
                readid = line.split()[0].split("_")[0]
                if readid not in bl_dic:
                    bl_dic[readid] = " "
            bamfile = pysam.AlignmentFile(f, "rb")
            out_bam = pysam.AlignmentFile(output_id, "wb", template = bamfile)
            for read in bamfile:
                readid = read.query_name.split("_")[0]
                if readid in bl_dic:
                    pass
                else: 
                    out_bam.write(read)
            bamfile.close()
            out_bam.close()
            os.system ("samtools index " + str(f)[0:-4] + "_bl.bam")  
            os.system ("{samtools} index {output_id}".format(
                samtools =  opt.samtools,
                output_id = output_id
            ))
