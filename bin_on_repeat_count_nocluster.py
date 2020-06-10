#! /usr/bin/env python
import sys
import os
import re
import subprocess
from optparse import OptionParser
from optparse import OptionGroup
import settings

if __name__ == "__main__":
    parser = OptionParser();
    group = OptionGroup(parser, "Main options")
    group.add_option("--im", dest = "m5dir", metavar = "[PATH]", help = "full path to  m5 folder [default = ./]")
    group.add_option("--ib", dest = "bamdir", metavar = "[PATH]", help = "full path to  bam folder [default = ./]")
    group.add_option("-o", dest = "outdir", metavar  ="[PATH]", help = "full path to output folder [default = ./")
    group.add_option("-p", default = settings.pbdagcon, dest = "pbdagcon", metavar = "[PATH]", help = "full path to pbgadcon binary [default pbdagcon in settings.py]")
    group.add_option("--param", default = settings.PBDAGCON_PARAM, dest = "param", metavar = "[STRING]", help = "pbdagcon parameters [default PBDAGCON_PARAM in settings.py]")
    group.add_option("--sa", default = settings.sambamba, dest = "sambamba", metavar = "[PATH]", help = "full path to sambamba binary [default sambamba in settings.py]")
    group.add_option("-b", default = settings.bwa, dest = "bwa", metavar = "[PATH]", help = "full path to bwa binary [default bwa in settings.py]")
    group.add_option("--rf", default = settings.full_ref, dest = "refgenome_full", metavar = "[PATH]", help = "full path to complete reference genome [default full_ref in settings.py]")
    group.add_option("-a", default = settings.INSERT, dest = "insert", metavar = "[STRING]", help = "name of insert chromosome [default INSERT in settings.py]")
    group.add_option("--cl", default = settings.MIN_CONS_LEN, dest = "cons_len", metavar = "INT", help = "minimum length (bp) for consensus calling [default MIN_CONS_LEN in settings.py]")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    pbdagcon_param = " -m " + str(opt.cons_len) + str(opt.param)

    if opt.outdir:
        outfolder = opt.outdir
    else:
        outfolder = os.getcwd()

    if not os.path.isdir("{outfolder}/tmp".format(outfolder=outfolder)):
        os.system("mkdir {outfolder}/tmp".format(outfolder=outfolder))
    if not os.path.isdir("{outfolder}/SH".format(outfolder=outfolder)):
        os.system("mkdir {outfolder}/SH".format(outfolder=outfolder))
    if not os.path.isdir("{outfolder}/bin_consensus".format(outfolder=outfolder)):
        os.system("mkdir {outfolder}/bin_consensus".format(outfolder=outfolder))
    if not os.path.isdir("{outfolder}/bin_consensus_folder".format(outfolder=outfolder)):
        os.system("mkdir {outfolder}/bin_consensus_folder".format(outfolder=outfolder))

    """Log GIT version + commit of repository"""
    if str(sys.argv[0]) == "python":
        repo = "/".join(sys.argv[1].split("/")[0:-1])
    else:
        repo = "/".join(sys.argv[0].split("/")[0:-1])
    os.system("git --git-dir={repo}/.git describe --tags > {outfolder}/GIT_bin_on_repeat_count.log".format(
        repo=repo,
        outfolder=outfolder
    ))
    os.system("git --git-dir={repo}/.git log --tags > {outfolder}/GIT_bin_on_repeat_count.log".format(
        repo=repo,
        outfolder=outfolder
    ))

    folders=os.listdir(opt.bamdir)
    countfolder=0
    for folder in folders:
        folder_id = str(folder).split("/")[-1]
        bamtarfile = "{bamdir}/{folder}/{folder_id}.tar".format(
            bamdir=opt.bamdir,
            folder=folder,
            folder_id=folder_id
        )
        m5tarfile = "{m5dir}/{folder}/{folder_id}.tar".format(
            m5dir=opt.m5dir,
            folder=folder,
            folder_id=folder_id
        )

        if os.path.isfile(bamtarfile):
            files=subprocess.getoutput("tar -tf {}".format(bamtarfile)).split()
            output_id = str(outfolder) + "/SH/Make_fasta_{}".format(countfolder)
            output_file = "{output_id}.sh".format(output_id=output_id) 
 
            if os.path.isfile(output_file): #check if file is already present. If so, do not continue.
                countfolder += 1
                continue

            for f in files:
                insert_count=0
                if "bai" not in str(f):
                    os.system("tar -axf {bamtarfile} {f}*".format(bamtarfile=bamtarfile, f=f))
                    insert_count=subprocess.getoutput("{sambamba} depth base -L {insert} {f} | cut -f3| sort -nk1 | tail -n1".format(
                        sambamba=opt.sambamba,
                        insert=opt.insert,
                        f=f
                    )).split("\n") #Last value is highest coverage.

                    if len(insert_count)> 1 and "failed" not in str(insert_count):
                        if "COV" not in insert_count[1]:
                            insert_count=insert_count[1]
                        else:
                            insert_count=0
                    else:
                        insert_count=0  
                    os.system("rm {}*".format(f))

                    if int(insert_count) >= 40:
                        action = "tar -axf {m5tarfile} {chopf}* -O | {pbdagcon} - {param} | sed 's/>/>{chopf}_/g' 1>> {outfolder}/bin_consensus_folder/{folder}_consensus_40.fasta".format(
                            m5tarfile=m5tarfile,
                            chopf=f[0:-11],
                            pbdagcon=opt.pbdagcon,
                            param=pbdagcon_param,
                            outfolder=outfolder,
                            folder=folder
                        )
                        os.system(action)

                    elif int(insert_count) > 0 and int(insert_count) <40:
                        for x in range(1, 40):
                            if int(insert_count) == x:
                                action = "tar -axf {m5tarfile} {chopf}* -O | {pbdagcon} - {param} | sed 's/>/>{chopf}_/g' 1>> {outfolder}/bin_consensus_folder/{folder}_consensus_{x}.fasta".format(
                                    m5tarfile=m5tarfile,
                                    chopf=f[0:-11],
                                    pbdagcon=opt.pbdagcon,
                                    param=pbdagcon_param,
                                    outfolder=outfolder,
                                    folder=folder,
                                    x=x 
                                )
                                os.system(action)
                    else: 
                        pass
            countfolder += 1  

    """ Merge fasta files """
    for x in range(1, 41):
         action = "find {outfolder}/bin_consensus_folder/ -iname \"*_consensus_{x}.fasta\" -exec cat {{}} \; >> {outfolder}/bin_consensus/consensus_{x}.fasta".format(
             outfolder=outfolder,
             x=x
         )
         os.system(action)

    """ Perform whole genome mapping """
    for x in range(1, 41):
        runid = "consensus_{x}".format(x=x)
        map_folder = "{outfolder}/bin_consensus/{runid}".format(outfolder=outfolder,runid=runid)
        action = "{bwa_setting} \"@RG\\tID:{runid}\\tSM:{runid}\\tPL:NANOPORE\\tLB:{runid}\" {refgenome_full} {map_folder}.fasta > {map_folder}_full_consensus.sam".format(
            bwa_setting=settings.BWA_MEM,
            runid=runid,
            refgenome_full=opt.refgenome_full,
            map_folder=map_folder
        )
        os.system(action)
 
        action = "{sambamba} view -S -f bam {map_folder}_full_consensus.sam > {map_folder}_full_consensus.bam\n".format(
            sambamba=opt.sambamba,
            map_folder=map_folder
        )
        os.system(action)

        action = "{sambamba} sort --tmpdir=./tmp {map_folder}_full_consensus.bam -o {map_folder}_full_consensus.sorted.bam\n".format(
            sambamba=opt.sambamba,
            map_folder=map_folder
        )
        os.system(action)

        os.system("rm {map_folder}_full_consensus.sam".format(map_folder=map_folder))
        os.system("rm {map_folder}_full_consensus.bam".format(map_folder=map_folder))

    """ Calculate allele count """
    os.chdir("{outfolder}/bin_consensus/".format(
        outfolder=outfolder
    ))   
 
    action = "source {venv} && {calculate} ".format(
        venv=settings.venv,
        calculate=settings.calculate
    )
    os.system(action)

    action = "rm {outfolder}/bin_consensus_folder/ -r".format(outfolder=outfolder)
    os.system(action)

