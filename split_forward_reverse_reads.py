#! /usr/bin/env python
import sys, os,re, commands
import pysam
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest="wkdir", metavar="[PATH]", help="full path to BAM input folder [default = ./")
        group.add_option("-o", dest="outdir", metavar="[PATH]", help="full path to output folder [default = ./")
        group.add_option("-b", default="/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa", dest="bwa", metavar="[PATH]", help="full path to bwa binary [default = /hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa ]")
        group.add_option("--sa", default="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba", dest="sambamba", metavar="[PATH]", help="full path to sambamba binary [default = /hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba]")
        group.add_option("--b5", default="/hpc/cog_bioinf/ridder/tools/bam2m5/bam2m5.py", dest="bam2m5", metavar="[PATH]", help="full path to bam2m5 binary [default = /hpc/cog_bioinf/ridder/tools/bam2m5/bam2m5.py]")
        group.add_option("--rf", default="/hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta", dest="refgenome_full", metavar="[PATH]", help="full path to complete reference genome [default = /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome//version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()


infolder=opt.wkdir # should be folder with sorted.bam
outfolder= opt.outdir

if outfolder.endswith('/'): # chop last "/" if present
    outfolder=outfolder[0:-1]

l=["forward","reverse"]
trim=0
bwa=str(opt.bwa) + " mem -t 2 -c 100 -M -R"
refgenome_full=opt.refgenome_full
sambamba=opt.sambamba
sample=outfolder.split("/")[-2].split("_")[0]

os.system ("mkdir "+str(outfolder)+"/tmp")
os.system("rm "+str(outfolder)+"/"+str(outfolder.split("/")[-1])+"_reverse.fasta")
os.system("rm "+str(outfolder)+"/"+str(outfolder.split("/")[-1])+"_forward.fasta")

for folder in os.listdir(infolder):
    tarfile = str(infolder)+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"
    if  os.path.isfile(tarfile):
        files=commands.getoutput("tar -tf "+str(tarfile)).split()
    for f in files:
        if "bai" not in f:
            if os.path.isfile(tarfile)== False:
                break
            os.system("tar -axf "+str(tarfile)+ " "+ str(f))
            in_file=pysam.AlignmentFile(str(outfolder)+"/"+str(f), "rb")
            forward=0
            reverse=0
            for line in in_file:
                if "TP53" in line.reference_name:
                    if line.flag == 16:
                        reverse+=1
                    else:
                        forward+=1
            os.system("rm "+ str(outfolder)+"/"+str(f))
            conf= f[0:-11]+str(".consensus")
            tarfilec = str(infolder.replace("bam", "consensus"))+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"
            if reverse>forward:
                os.system("tar -axf "+str(tarfilec)+ " "+ str(conf) + " -O >> "+ str(outfolder.split("/")[-1])+"_reverse.fasta")
            elif forward>reverse:
                os.system("tar -axf "+str(tarfilec)+ " "+ str(conf) + " -O >> "+ str(outfolder.split("/")[-1])+"_forward.fasta")
            else:
                pass # discard all reads that are 50/50 f/r      


for item in l:
    os.system (str(bwa) + "\"@RG\\tID:"+str(item)+"\\tSM:"+(sample)+"\\tPL:NANOPORE\\tLB:"+str(sample)+"\" "+str(refgenome_full)+" "+str(outfolder.split("/")[-1])+"_"+str(item)+".fasta | "+str(sambamba)+ " view -S -f bam /dev/stdin | "+ str(sambamba)+ " sort -t 2 --tmpdir=./tmp /dev/stdin -o "+str(outfolder.split("/")[-1])+"_"+str(item)+ ".sorted.bam")
    os.system("rm "+str(outfolder.split("/")[-1])+"_"+str(item)+".fasta")

