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
    group.add_option("--id", dest = "sampleid", metavar = "[STRING]", help = "sampleID")
    group.add_option("--mp", default =settings.SPLIT_MIN_PERC, dest = "minperc", metavar = "[FLOAT]", help = "minimum percentage of forward or reverse reads required to include in forward or reverse bin [default SPLIT_MIN_PERC in settings.py")
    group.add_option("--sa", default = settings.sambamba, dest = "sambamba", metavar = "[PATH]", help = "full path to sambamba executable [default sambamba in settings.py]")
    group.add_option("--b5", default = settings.bam2m5, dest = "bam2m5", metavar = "[PATH]", help = "full path to bam2m5 executable [default bam2m5 in settings.py")
    group.add_option("--rf", default = settings.full_ref, dest = "refgenome_full", metavar = "[PATH]", help = "full path to complete reference genome [default full_ref in settings.py]")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    trim = settings.TRIM
    states = settings.STATES

    if opt.outdir.endswith('/'): # chop last "/" if present
        outfolder = opt.outdir[0:-1]
    else:
        outfolder = opt.outdir

    if opt.sampleid:
        sample = opt.sampleid
    else:
        sample = outfolder.split("/")[-2].split("_")[0]

    run = str(outfolder.split("/")[-1])

    if opt.gene:
        gene = str(opt.gene)
    else:
        gene = "TP53"

    # Log GIT version + commit of repository
    if str(sys.argv[0]) == "python":
        repo = "/".join(sys.argv[1].split("/")[0:-1])
    else:
        repo = "/".join(sys.argv[0].split("/")[0:-1])
    os.system("git --git-dir={repo}/.git describe --tags > {outfolder}/GIT_split_forward_reverse_reads.log".format(
        repo=repo,
        outfolder=outfolder
    ))

    os.system("git --git-dir={repo}/.git log --tags > {outfolder}/GIT_split_forward_reverse_reads.log".format(
        repo=repo,
        outfolder=outfolder
    ))

    if not os.path.isdir("{outfolder}/tmp".format(outfolder=outfolder)):
        os.system("mkdir {outfolder}/tmp".format(outfolder=outfolder))

    forward_fasta = "{outfolder}/{run}_forward.fasta".format(outfolder=outfolder, run=run)
    reverse_fasta = "{outfolder}/{run}_reverse.fasta".format(outfolder=outfolder, run=run)

    if os.path.isfile(forward_fasta):
        os.system("rm {forward_fasta}".format(forward_fasta=forward_fasta))
    if os.path.isfile(reverse_fasta):
        os.system("rm {reverse_fasta}".format(reverse_fasta=reverse_fasta))

    for folder in os.listdir(opt.wkdir):
        folder_id = str(folder).split("/")[-1]
        tarfile = "{wkdir}/{folder}/{folder_id}.tar".format(
            wkdir=opt.wkdir,
            folder=folder,
            folder_id=folder_id
        )

        if os.path.isfile(tarfile): # skip empty folders
            files = subprocess.getoutput("tar -tf {}".format(tarfile)).split()
            for f in files:
                if "bai" not in f:
                    f_path = "{outfolder}/{f}".format(outfolder=outfolder,f=f)
                    if os.path.isfile(tarfile) == False:
                        break
                    os.system("tar -xf {tarfile} {f}".format(
                        tarfile=tarfile,
                        f=f
                    ))
                    in_file = pysam.AlignmentFile(f_path, "rb")
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
                    os.system("rm {}".format(f_path))
                    conf = "{}.consensus".format(f[0:-11])
                    tarfilec = "{cons_wkdir}/{folder}/{folder_id}.tar".format(
                        cons_wkdir=opt.wkdir.replace("bam", "consensus"),
                        folder=folder,
                        folder_id=folder_id
                    )
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
                        if rev_readperc > float(opt.minperc):
                            os.system("tar -xOf {tarfilec} {conf} >> {reverse_fasta}".format(
                                tarfilec=tarfilec,
                                conf=conf,
                                reverse_fasta=reverse_fasta
                            ))

                            count_rev_new += 1
                        elif for_readperc > float(opt.minperc):
                            os.system("tar -xOf {tarfilec} {conf} >> {forward_fasta}".format(
                                tarfilec=tarfilec,
                                conf=conf,
                                forward_fasta=forward_fasta
                            ))
                            count_for_new += 1
                        else:
                            pass # discard reads
                    else:
                        rev_readperc = 0
                        for_readperc = 0
                    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(f, reverse, forward, count_rev_new, count_for_new, rev_readperc, for_readperc))

    for state in states:
        action = "{bwa} \"@RG\\tID:{state}\\tSM:{sample}\\tPL:NANOPORE\\tLB:{sample}\" {refgenome_full} {run}_{state}.fasta | {sambamba} view -S -f bam /dev/stdin | {sambamba} sort -t 2 --tmpdir=./tmp /dev/stdin -o {run}_{state}.sorted.bam".format(
            bwa=settings.BWA_MEM,
            state=state,
            sample=sample,
            refgenome_full=opt.refgenome_full,
            run=run,
            sambamba=opt.sambamba
        )
        os.system(action)
        os.system("{calculate}\n".format(calculate=settings.calculate))
        os.system("rm {run}_{state}.fasta".format(run=run, state=state))
        os.system("rm -r tmp".format(run=run, state=state))
