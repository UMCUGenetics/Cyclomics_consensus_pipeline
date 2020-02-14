#! /usr/bin/env python
import sys, os, re
import glob
import commands
import gzip

from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest="wkdir", metavar="[PATH]", help="full path to raw data folder [default = ./")
        group.add_option("-o", dest="outdir", metavar="[PATH]", help="full path to output folder [default = ./")
        group.add_option("-s", default=4000, dest="size", metavar="[INT]", help="number of reads per file [default = 4000]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

if opt.wkdir:
    wkdir=opt.wkdir
else:
   sys.exit("provide folder to raw data" )

if opt.outdir:
    outdir=opt.outdir
else:
    sys.exit("provide output folder" )

size=int(opt.size)

files= commands.getoutput("find "+wkdir+ " -maxdepth 1 -iname \"*fastq\"").split()

if len(files)==0: ## likely that fastq are compressed
    files= commands.getoutput("find "+wkdir+ " -maxdepth 1 -iname \"*fastq.gz\"").split()
    gz="on"
else:
    gz="off"

os.system("mkdir "+str(outdir))

for f in files:
    if gz=="off":
        lines=open(f, "r").readlines()
        file_name=f.split("/")[-1][0:-6]
    else:
        lines=gzip.open(f, "rb").readlines()
        file_name=f.split("/")[-1][0:-9]

    print file_name
    x=0
    c=0
    bin = 0
    p_list=[]
    while x< len(lines):
        if "@" in lines[x]:     # bypass nasty bug in which fastq's are corrupted. Corrupted will be skipped, not fixed
            splitline=lines[x].split()
            #read=str(splitline[0].split("@")[1])+"_read_"+splitline[3].split("=")[1]+"_ch_"+splitline[4].split("=")[1]
            p_list+=[[lines[x],lines[x+1],lines[x+2],lines[x+3]]]
            x+=4
            c+=1
        else:
            x+=1          
        if c == size:
            write_file=open(str(outdir)+"/"+str(file_name)+"_"+str(bin)+".fastq","w")
            for item in p_list:
                for readline in item:
                    write_file.write(readline.rstrip()+"\n")
            write_file.close()
            os.system("bgzip "+str(outdir)+"/"+str(file_name)+"_"+str(bin)+".fastq")
            p_list=[]
            c=0
            bin+=1

    write_file=open(str(outdir)+"/"+str(file_name)+"_"+str(bin)+".fastq","w") 
    for item in p_list:
        for readline in item:
            write_file.write(readline.rstrip()+"\n")
    write_file.close()
    os.system("bgzip "+str(outdir)+"/"+str(file_name)+"_"+str(bin)+".fastq")


