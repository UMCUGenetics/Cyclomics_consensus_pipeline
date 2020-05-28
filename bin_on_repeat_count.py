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
    group.add_option("--project", default = settings.project, dest = "project", metavar = "STRING", help = "SGE project for submitting jobs [default project in settings.py]")
    group.add_option("-m", default = settings.mail, dest = "mail", metavar = "[STRING]", help = "email used for job submitting [default mail in settings.py]")
    group.add_option("-a", default = settings.INSERT, dest = "insert", metavar = "[STRING]", help = "name of insert chromosome [default INSERT in settings.py]")
    group.add_option("--tlow", default = settings.SLURM_JOB_TIME_LOW, dest = "timeslotlow", metavar = "[TIME]", help = "time slot for jobs [default SLURM_JOB_TIME in settings.py]")
    group.add_option("--tmed", default = settings.SLURM_JOB_TIME_MED, dest = "timeslotmed", metavar = "[TIME]", help = "time slot for jobs [default SLURM_JOB_TIME in settings.py]")
    group.add_option("--maxfiles", dest = "maxfiles", metavar = "[INT]", help = "maximum number of files to be processed [default = all]")
    group.add_option("--mem_full", default = settings.MAX_MEM_FULL, dest = "max_mem_full", metavar = "[INT]", help = "memory used for jobs [default MAX_MEM_FULL in settings.py]")
    group.add_option("--mem_target", default = settings.MAX_MEM_TARGET, dest = "max_mem_target", metavar = "[INT]", help = "memory used for jobs [default MAX_MEM_TARGET in settings.py]")
    group.add_option("--threads", default = settings.THREADS, dest = "threads", metavar = "[INT]", help = "number threads used for jobs [default THREADS in settings.py]")
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

    test=[]
    job_id=[]
    folders=os.listdir(opt.bamdir)
    if opt.maxfiles:
        maxfiles=int(opt.maxfiles)
    else:
        maxfiles=int(len(folders))

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

        if countfolder >= maxfiles:
            break 
        if os.path.isfile(bamtarfile):
            files=subprocess.getoutput("tar -tf {}".format(bamtarfile)).split()
            output_id = str(outfolder) + "/SH/Make_fasta_{}".format(countfolder)
            output_file = "{output_id}.sh".format(output_id=output_id) 
 
            if os.path.isfile(output_file):
                test += [output_file]
                countfolder += 1
                continue
            write_file=open(output_file,"w")
            write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project}\n#SBATCH --mem={mem}G\n#SBATCH --export=NONE\n#SBATCH -o {output_id}.output -e {output_id}.error \n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={mail}\n".format(
                timeslot=opt.timeslotlow,
                project=opt.project,
                mem=opt.max_mem_target,
                output_id=output_id,
                mail=opt.mail
            ))

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
                        write_file.write("tar -axf {m5tarfile} {chopf1}* -O | {pbdagcon} - {param} | sed 's/>/>{chopf2}_/g' 1> {outfolder}/bin_consensus_folder/{folder}_consensus_40.fasta\n".format(
                            m5tarfile=m5tarfile,
                            chopf1=f[0:-11],
                            pbdagcon=opt.pbdagcon,
                            param=pbdagcon_param,
                            chopf2=f[0:-10],
                            outfolder=outfolder,
                            folder=folder
                        ))
                    elif int(insert_count) > 0 and int(insert_count) <40:
                        for x in range(1, 40):
                            if int(insert_count) == x:
                                write_file.write("tar -axf {m5tarfile} {chopf1}* -O | {pbdagcon} - {param} | sed 's/>/>{chopf2}_/g' 1> {outfolder}/bin_consensus_folder/{folder}_consensus_{x}.fasta\n".format(
                                    m5tarfile=m5tarfile,
                                    chopf1=f[0:-11],
                                    pbdagcon=opt.pbdagcon,
                                    param=pbdagcon_param,
                                    chopf2=f[0:-10],
                                    outfolder=outfolder,
                                    folder=folder,
                                    x=x 
                                 ))
                    else: 
                        pass

            write_file.close()
            test += [output_file] 
        countfolder += 1  

    folder_dic={}
    for item in test:
        folder="_".join(item.split("/")[-1].split("_")[0:-1])
        if folder not in folder_dic:
            folder_dic[folder] = 1
        else:
            folder_dic[folder] += 1

    for item in folder_dic:
        out_file_id = "{outfolder}/SH/{item}".format(outfolder=outfolder,item=item)
        out_file = "{out_file_id}_array.sh".format(out_file_id=out_file_id)
        array_folder_id = "{outfolder}/SH/{folder}".format(outfolder=outfolder,folder=folder)

        write_file=open(out_file,"w")
        write_file.write("#!/bin/bash\n#SBATCH -t {timeslot}\n#SBATCH --mem={mem}G\n#SBATCH --account={project}\n#SBATCH --mail-type=FAIL\n#SBATCH --export=NONE\n#SBATCH --mail-user={mail}\n#SBATCH -o {array_folder_id}_%A_%a.output\n#SBATCH -e {array_folder_id}_%A_%a.error\n#SBATCH --array=0-{number}%4\n".format(
                timeslot=opt.timeslotmed,
                mem=opt.max_mem_target,
                project=opt.project,
                mail=opt.mail,
                output_id=output_id,
                array_folder_id=array_folder_id,
                number=int(folder_dic[folder])-1
            ))

        write_file.write("sh {out_file_id}_$SLURM_ARRAY_TASK_ID\.sh\n".format(out_file_id=out_file_id))
        write_file.close()
        job_output = subprocess.getoutput("sbatch {}".format(out_file))
        job_id_make_fasta = job_output.split()[3]


    merge_file = "{outfolder}/SH/merge_fasta.sh".format(outfolder=outfolder)
    write_file=open(merge_file,"w")
    write_file.write("#!/bin/bash\n#SBATCH -t {timeslot}\n#SBATCH --export=NONE\n#SBATCH --account={project}\n#SBATCH --mem={mem}G\n#SBATCH -o {merge_file}.output -e {merge_file}.error \n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={mail}\n".format(
        timeslot=opt.timeslotlow,
        project=opt.project,
        mem=opt.max_mem_target,
        merge_file=merge_file,
        mail=opt.mail
    ))

    for x in range(1, 41):
         write_file.write("find {outfolder}/bin_consensus_folder/ -iname \"*_consensus_{x}.fasta\" -exec cat {{}} \; >> {outfolder}/bin_consensus/consensus_{x}.fasta\n".format(
             outfolder=outfolder,
             x=x
         ))
    write_file.close()

    if len(job_id) >500:	## this will set the maximum hold jobs to 500 (max HPC slurm is 1000)
        job_id=job_id[-500:]
    job_output=subprocess.getoutput("sbatch -c 2 --depend={job_id_make_fasta} {merge_file}".format(job_id_make_fasta=job_id_make_fasta, merge_file=merge_file)) 
    job_id_merge=job_output.split()[3]

    test=[]
    def write_new_file(x,test):
        runid = "consensus_{x}".format(x=x)
        map_file = "{outfolder}/SH/{runid}_mapping.sh".format(outfolder=outfolder,runid=runid)
        write_file=open(map_file,"w")
        write_file.write("#!/bin/bash\n#SBATCH -t {timeslot}\n#SBATCH --export=NONE\n#SBATCH --account={project}\n#SBATCH --mem={mem}G\n#SBATCH -o {map_file}.output -e {map_file}.error \n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={mail}\n".format(
            timeslot=opt.timeslotlow,
            project=opt.project,
            mem=opt.max_mem_target,
            map_file=map_file,
            mail=opt.mail
        ))

        map_folder = "{outfolder}/bin_consensus/{runid}".format(outfolder=outfolder,runid=runid)
        write_file.write("{bwa_setting} \"@RG\\tID:{runid}\\tSM:{runid}\\tPL:NANOPORE\\tLB:{runid}\"{refgenome_full} {map_folder}.fasta > {map_folder}_full_consensus.sam\n".format(
            bwa_setting=settings.BWA_MEM,
            runid=runid,
            refgenome_full=opt.refgenome_full,
            map_folder=map_folder
        ))


        write_file.write("{sambamba} view -S -f bam {map_folder}_full_consensus.sam > {map_folder}_full_consensus.bam\n".format(
            sambamba=opt.sambamba,
            map_folder=map_folder
        )) 

        write_file.write("{sambamba} sort -t {threads} --tmpdir=./tmp {map_folder}_full_consensus.bam -o {map_folder}_full_consensus.sorted.bam\n".format(
            sambamba=opt.sambamba,
            threads=opt.threads,
            map_folder=map_folder
        ))
        write_file.write("sleep 2\n")
        write_file.write("rm {map_folder}_full_consensus.sam\nsleep 2\n)".format(map_folder=map_folder))
        write_file.write("rm {map_folder}_full_consensus.bam\nsleep 2\n)".format(map_folder=map_folder))
        write_file.write("mv {map_folder}_full_consensus.sorted.bam* {outfolder}/bin_consensus/\nsleep 2\n".format(map_folder=map_folder, outfolder=outfolder))
        write_file.close()
        test += [map_file]
        return test

    for x in range(1, 41):
        job_id = write_new_file(x, job_id)
        test = write_new_file(x, test)

    cons_file = "{outfolder}/SH/consensus_calling_array.sh".format(outfolder=outfolder) 
    write_file = open(cons_file,"w")
    write_file.write("#!/bin/bash\n#SBATCH -t {timeslot}\n#SBATCH --mem={mem}G\n#SBATCH --account={project}\n#SBATCH --mail-type=FAIL\n#SBATCH --export=NONE\n#SBATCH --mail-user={mail}\n \
                     #SBATCH -o {cons_file}_%A_%a.output\n#SBATCH -e {cons_file}_%A_%a.error\n#SBATCH --array=0-{number}%8\n".format(
        timeslot=opt.timeslotmed,
        mem=opt.max_mem_target,
        project=opt.project,
        mail=opt.mail,
        cons_file=cons_file,
        number=len(test)
    )) 

    write_file.write("sh {outfolder}/SH/consensus__$SLURM_ARRAY_TASK_ID\_mapping.sh\n".format(outfolder=outfolder)) 
    write_file.close()
    job_output=subprocess.getoutput("sbatch --depend={job_id_merge} {cons_file}".format(job_id_merge=job_id_merge, cons_file=cons_file))
