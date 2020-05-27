#! /usr/bin/env python
import sys, os, re
import glob
import subprocess
import time
import gzip
from optparse import OptionParser
from optparse import OptionGroup
import settings


if __name__ == "__main__":
    parser = OptionParser();
    group = OptionGroup(parser, "Main options")
    group.add_option("-i", dest = "wkdir", metavar = "[PATH]", help = "full path to FASTQ folder [default = ./")
    group.add_option("-o", dest = "outdir", metavar = "[PATH]", help = "full path to output folder [default = ./")
    group.add_option("-b", dest = "blacklist", metavar = "[PATH]", help = "full path to blacklist file [default = off]")
    group.add_option("-c", default = settings.DAGCON_MIN_COV, dest = "coverage", metavar = "[INT]", help = "minimum coverage required for assembly [default = pgdagcon_coverage in settings.py]")
    group.add_option("-m", default = settings.mail, dest = "mail", metavar = "[STRING]", help = "email used for job submitting [default mail in settings.py]")
    group.add_option("-t", default = settings.SLURM_JOB_TIME_LOW, dest = "timeslot", metavar = "[TIME]", help = "time slot for jobs [deafult SLURM_JOBTIME in settings.py]")
    group.add_option("-n", default = settings.MAX_READS_JOB, dest = "number", metavar="[INT]", help = "max number of reads within a scatterjob [default MAX_READS_JOB in settings.py]")
    group.add_option("-p", dest = "prefix", metavar = "[STRING]", help = "prefix of BAM names [default = off (FASTQ input folder name)]")
    group.add_option("-s", default = settings.SLURM_PARALLEL_JOBS, dest = "slurm_job", metavar="[INT]", help = "number of parallel slurm [default SLURM_PARALLEL_JOPS in settings.py]")
    group.add_option("--mem_target", default = settings.MAX_MEM_TARGET, dest = "max_mem_target", metavar = "[INT]", help = "memory used for jobs [default MAX_MEM_TARGET in setttings.py]")
    group.add_option("--mem_full", default = settings.MAX_MEM_FULL, dest = "max_mem_full", metavar = "[INT]", help = "memory used for jobs [default MAX_MEM_FULL in settings.py]")
    group.add_option("--threads", default = settings.THREADS, dest = "threads", metavar = "[INT]", help = "number threads used for jobs [default THREADS in settings.py]")
    group.add_option("--cl", default = settings.MIN_CONS_LEN, dest = "cons_len", metavar = "INT", help = "minimum length (bp) for consensus calling [default MIN_CONS_LEN in settings.py]")      
    group.add_option("--project", default = settings.project, dest = "project", metavar = "STRING", help = "SGE project for submitting jobs [default project in settings.py]")
    group.add_option("--bwa", default = settings.bwa, dest = "bwa", metavar = "[PATH]", help = "full path to bwa executable [default bwa in settings.py ]")
    group.add_option("--sa", default = settings.sambamba, dest = "sambamba", metavar = "[PATH]", help = "full path to sambamba executable [default sambamba in settings.py]")

    group.add_option("--la", default = settings.lastal, dest = "lastal", metavar = "[PATH]", help = "full path to lastal executable [default last in settings.py]")
    group.add_option("--ls", default = settings.lastsplit, dest = "lastsplit", metavar = "[PATH]", help = "full path to last-split executable [default last in settings.py]")
    group.add_option("--lp", default = settings.lastparam, dest = "lastparam", metavar = "[PATH]", help = "full path to last-param file [default last in settings.py]")
    group.add_option("--maf", default = settings.mafconvert, dest = "mafconvert", metavar = "[PATH]", help = "full path to maf-convert executable [default last in settings.py]")

    group.add_option("-e", default = settings.venv, dest = "venv", metavar = "[ENV]", help = "full path to python enviroment [default venv in settings.py]")
    group.add_option("--b5", default = settings.bam2m5, dest = "bam2m5", metavar = "[PATH]", help = "full path to bam2m5 executable [default bam2m5 in settings.py]")

    group.add_option("--pbdagcon", default = settings.pbdagcon, dest = "pbdagcon", metavar = "[PATH]", help = "full path to pbgadcon executable [default pbdagcon in settings.py]")
    group.add_option("--rf", default = settings.full_ref, dest = "refgenome_full", metavar = "[PATH]", help = "full path to complete reference genome [default full_ref in settings.py]")
    group.add_option("--rt", default = settings.target_ref, dest = "refgenome_target", metavar = "[PATH]", help = "full path to targeted reference genome [default = target_ref in settings.py]")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    trim = settings.TRIM

    if not opt.wkdir:
        sys.exit("provide input folder")
    wkdir = opt.wkdir

    if opt.outdir:
        outdir = opt.outdir
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    else:
        sys.exit("provide output folder")

    if opt.prefix:
        runid = str(opt.prefix)
    else:
        runid = wkdir.split("/")[-2]

    if "fasta" not in opt.refgenome_target:
        sys.exit("please provide fasta file for target reference genome")
    if "fasta" not in opt.refgenome_full:
        sys.exit("please provide fasta file for full reference genome")

    dic_bl = {}
    if opt.blacklist:
        blacklist = open(opt.blacklist,"r").readlines()
        for item in blacklist:
            if "read" in item:
                if item.split()[1].split("_")[1] not in dic_bl:
                    dic_bl[item.split()[1].split("_")[1]] = "_".join(item.split()[2:])
                else:
                    #print("Warning, read found multiple times in blacklist {}".format(item.split()[1].split("_")[1]))
                    pass

    """Log GIT version + commit of repository"""
    if str(sys.argv[0]) == "python":
        repo = "/".join(sys.argv[1].split("/")[0:-1])
    else:
        repo = "/".join(sys.argv[0].split("/")[0:-1])

    os.system("git --git-dir={repo}/.git describe --tags > {outdir}/GIT_run_dagcon_consensus.log".format(repo = repo, outdir = outdir))
    os.system("git --git-dir={repo}/.git log >> {outdir}/GIT_run_dagcon_consensus.log".format(repo = repo, outdir = outdir))

    """Print log of used parameters"""
    write_file = open(outdir + "/" + str(runid) + ".log", "w")
    (options, args) = parser.parse_args()
    for item in vars(options):
        write_file.write("{item}\t{varsitem}\n".format(item = item, varsitem = vars(options)[item]))
    write_file.close()

    if wkdir.endswith('/'): # chop last "/" if present
        wkdir = wkdir[0:-1]

    if os.path.isdir(outdir + "/bam") or os.path.isdir(outdir + "/m5") or os.path.isdir(outdir + "/consensus"):
        sys.exit("Please remove previous analysis first")

    refgenome_target_db = opt.refgenome_target[0:-6]
    refgenome_target = opt.refgenome_target
    files = subprocess.getoutput("find {} -iname \"*fastq\"".format(wkdir)).split()
    if len(files) == 0: ## likely that fastq are compressed
        files = subprocess.getoutput("find {} -iname \"*fastq.gz\"".format(wkdir)).split()
        gz = "on"    
    else:
        gz = "off"
    fastq_dic = {}
    line_dic = {}
    job_list = []
    for fastq_files in files:
        if gz == "off":
             lines = open(fastq_files, "r").readlines()
        else:
            lines = gzip.open(fastq_files, "rb").readlines()
        fastqline = 0
        """Make dictionary of line-position of reads within the fastq file. This will be used to speed up targeted mapping"""
        while fastqline < len(lines):
            if "@".encode('UTF-8') in lines[fastqline]:	# bypass nasty bug in which fastq's are corrupted. Corrupted will be skipped, not fixed
                splitline = lines[fastqline].split()
                read_id = splitline[0].split("@".encode('UTF-8'))[1].decode("utf-8")
                read = splitline[2].split("=".encode('UTF-8'))[1].decode("utf-8")
                channel = splitline[3].split("=".encode('UTF-8'))[1].decode("utf-8")
                read_new = "{0}_read_{1}_ch_{2}".format(read_id, read, channel)
                if str(fastq_files) not in line_dic:
                    line_dic[str(fastq_files)] = [{read_new:[fastqline + 1, fastqline + 4]}]
                else:
                    line_dic[str(fastq_files)] += [{read_new:[fastqline + 1, fastqline + 4]}]
                fastqline += 4
            else: # loop until next @ is found 
                fastqline += 1

    total = 0
    subnumber = 0
    os.system("mkdir {}/SH".format(outdir))
    for fastq_file in line_dic:
        if gz == "off":  #fastq not gz
            if "pass" in fastq_file or "fail" in fastq_file:  # handles run with or without pass/fail reads
                folder = "{0}_{1}".format(fastq_file.split("/")[-2], fastq_file.split("/")[-1][0:-6])
            else:
                folder = fastq_file.split("/")[-1][0:-6]
        elif gz == "on":  #fastq.gz
            if "pass" in fastq_file or "fail" in fastq_file: # handles run with or without pass/fail reads
                folder = "{0}_{1}".format(fastq_file.split("/")[-2], fastq_file.split("/")[-1][0:-9])
            else:
                folder = fastq_file.split("/")[-1][0:-9]
        count = 0

        """Make targeted mapping jobs for each fastq file, split by MAX_READS_JOB"""
        write_file = open("{0}/SH/Targetmapping_{1}.sh".format(outdir, subnumber), "w")
        write_file.write("#!/bin/bash\n#SBATCH -t {time} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {outdir}/SH/{folder}_job_{subnumber}.output\n#SBATCH -e {outdir}/SH/{folder}_job_{subnumber}.error \n#SBATCH --mail-user={mail}\n".format(
            time=opt.timeslot ,
            project=opt.project,
            mem=opt.max_mem_target ,
            outdir=outdir ,
            folder=folder, 
            subnumber=subnumber,
            mail=opt.mail
        ))

        write_file.write(". {}\n".format(opt.venv))
        write_file.write("echo \"Start poststats    \" `date` \"    \" `uname -n`\n\n")
        os.system("mkdir -p {0}/bam/{1}_{2}".format(outdir, folder, subnumber))
        os.system("mkdir -p {0}/m5/{1}_{2}".format(outdir, folder, subnumber))
        os.system("mkdir -p {0}/consensus/{1}_{2}".format(outdir, folder, subnumber))
        for f in line_dic[fastq_file]:
            total += 1
            position = list(f.values())[0]
            fastq = list(f.keys())[0]
            if opt.blacklist:
                if fastq.split("_")[0] in dic_bl:
                    if "OK" not in dic_bl[fastq.split("_")[0]]:
                        readf = "fail"
                    else:
                        readf = "pass" 
                else:
                    readf = "pass"  
            else: # if no blacklist, all reads are passed
                 readf = "pass"

            if readf == "pass": 
                """Make BAM file."""
                if gz == "on":
                    cat = "zcat" 
                else:
                    cat = "cat"          

                write_file.write("{cat} {item} | sed -n \'{pos0},{pos1}\'p | {lastal} -Q 1 -p {lastparam} {refgenome_target_db} -P 1 /dev/stdin | {lastsplit} | {mafconvert} -f {refgenome_target_db}.dict sam -r \"ID:{fastq} PL:nanopore SM:{fastq}\" /dev/stdin | {sambamba} view -S -f bam /dev/stdin | {sambamba} sort -t 1 /dev/stdin -o {outdir}/bam/{folder}_{subnumber}/{fastq}.sorted.bam\n ".format(
                        cat=cat, 
                        item=fastq_file,
                        pos0=position[0],
                        pos1=position[1],
                        lastal=opt.lastal,
                        lastparam=opt.lastparam,
                        refgenome_target_db=refgenome_target_db,
                        lastsplit=opt.lastsplit,
                        mafconvert=opt.mafconvert,
                        fastq=fastq,
                        sambamba=opt.sambamba,
                        outdir=outdir,
                        folder=folder,
                        subnumber=subnumber,
                        ))          

                """Add Bam to TAR ball"""
                write_file.write("cd {outdir}/bam/{folder}_{subnumber}\n".format(
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber
                    ))
                write_file.write("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.sorted.bam*\n".format(
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq
                    ))

                """Extract BAM file from TAR ball and stdout to command to make m5"""
                write_file.write("tar -axf {outdir}/bam/{folder}_{subnumber}/{folder}_{subnumber}.tar {fastq}.sorted.bam -O | python {bam2m5} - {refgenome_target} {outdir}/m5/{folder}_{subnumber}/{fastq}.sorted.m5\n".format(
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq,
                    bam2m5=opt.bam2m5,
                    refgenome_target=refgenome_target
                    ))

                """Add m5 to TAR ball"""
                write_file.write("cd {outdir}/m5/{folder}_{subnumber}\n".format(
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber
                    ))
                write_file.write("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.sorted.m5*\n".format(
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq
                    ))

                """Extract M5 file from TAR ball and stdout to command"""
                write_file.write("echo {fastq} >> {outdir}/consensus/Full_insert_count_{folder}_{subnumber}.txt\n".format(
                    fastq=fastq,
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber
                    ))

                write_file.write("tar -axf {outdir}/m5/{folder}_{subnumber}/{folder}_{subnumber}.tar {fastq}.sorted.m5 -O |{pbdagcon} - -m {cons_len} -c {coverage} -t {trim} -j {threads} 1> {outdir}/consensus/{folder}_{subnumber}/{fastq}.consensus 2>> {outdir}/consensus/Full_insert_count_{folder}_{subnumber}.txt\n".format(
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq,
                    pbdagcon=opt.pbdagcon,
                    cons_len=opt.cons_len,
                    coverage=opt.coverage,
                    trim=trim,
                    threads=opt.threads
                    ))
                write_file.write("sed -i \'s/>/>{fkey}_/g\' {outdir}/consensus/{folder}_{subnumber}/{fastq}.consensus\n".format(
                    fkey=list(f.keys())[0],
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq
                    ))
 
                """Add consensus to TAR ball"""
                write_file.write("cd {outdir}/consensus/{folder}_{subnumber}\n".format(
                    outdir=outdir,
                    folder=folder,
                    subnumber=subnumber
                    ))
                write_file.write("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.consensus*\n".format(
                    folder=folder,
                    subnumber=subnumber,
                    fastq=fastq
                    ))

            if count == (opt.number):
                job_list += ["{0}/SH/Targetmapping_{1}.sh".format(outdir, subnumber)]
                write_file.close()
                subnumber += 1
                if total == len(line_dic[fastq_file]):
                    pass
                else:
                    write_file = open("{0}/SH/Targetmapping_{1}.sh".format(outdir, subnumber), "w")
                    write_file.write(". {}\n".format(opt.venv))
                    write_file.write("echo \"Start poststats    \" `date` \"    \" `uname -n`\n")
                    os.system("mkdir -p {0}/bam/{1}_{2}".format(outdir, folder, subnumber))
                    os.system("mkdir -p {0}/m5/{1}_{2}".format(outdir, folder, subnumber))
                    os.system("mkdir -p {0}/consensus/{1}_{2}".format(outdir, folder, subnumber))
                count = 0
            else:
                count += 1

        job_list += ["{0}/SH/Targetmapping_{1}.sh".format(outdir, subnumber)]
        subnumber += 1
        write_file.close()

    """Submit targeted mapping jobs as array to Slurm"""
    for item in job_list:
        write_file=open(str(outdir) + "/SH/Targetmapping_array.sh","w")
        write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \
            \n#SBATCH --mem={mem}G \
            \n#SBATCH -o {outdir}/SH/Targetmapping_%A_%a.output \
            \n#SBATCH -e {outdir}/SH/Targetmapping_%A_%a.error \
            \n#SBATCH --mail-type=FAIL \
            \n#SBATCH --account={project} \
            \n#SBATCH --export=NONE \
            \n#SBATCH --mail-user={mail} \
            \n#SBATCH --array=0-{files}%{slurm}\n".format(
            timeslot=opt.timeslot,
            mem=opt.max_mem_target,
            outdir=outdir,
            project=opt.project,
            mail=opt.mail,
            files=subnumber - 1,
            slurm=opt.slurm_job
        ))
        write_file.write("sh " + str(outdir) + "/SH/Targetmapping_$SLURM_ARRAY_TASK_ID\.sh\n")
        write_file.close()

    job_output=subprocess.getoutput("sbatch " + str(outdir) + "/SH/Targetmapping_array.sh")
    job_id_make_fasta = job_output.split()[3]

    """Make full genome mapping job and submit"""
    write_file = open("{}/SH/full_target_mapping.sh".format(outdir), "w")
    write_file.write("#!/bin/bash\n#SBATCH -t {time} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {outdir}/SH/full_target_mapping.output\n#SBATCH -e {outdir}/SH/full_target_mapping.error \n#SBATCH --mail-user={mail}\n".format(
            time=opt.timeslot ,
            project=opt.project,
            mem=opt.max_mem_full ,
            outdir=outdir ,
            mail=opt.mail
        ))
    write_file.write("mkdir {}/tmp\n".format(outdir))
    """Extract all file from all consensus TAR balls and put in FASTA file"""
    write_file.write("find {outdir}/consensus/ -type f -iname \"*tar\" -exec tar -O -xf {{}} \; >> {outdir}/{runid}_full_consensus.fasta\n".format(
        runid = runid,
        outdir = outdir
        ))  
    write_file.write("{bwa} mem -t 2 -c 100 -M -R \"@RG\\tID:{runid}\\tSM:{runid}\\tPL:NANOPORE\\tLB:{runid}\" {refgenome_full} {outdir}/{runid}_full_consensus.fasta | {sambamba} view -S -f bam /dev/stdin | {sambamba} sort -t 2 --tmpdir={outdir}/tmp /dev/stdin -o {outdir}/{runid}_full_consensus.sorted.bam\n".format(
        bwa=opt.bwa,
        runid=runid,
        refgenome_full=opt.refgenome_full,
        outdir=outdir,
        sambamba=opt.sambamba
        ))
    write_file.write("gzip {outdir}/{runid}_full_consensus.fasta\n".format(outdir = outdir, runid = runid))
    write_file.write("cd {}/SH\n".format(outdir))
    write_file.write("zip -m SH.zip *\n")
    write_file.close()

    #os.system("sbatch -c 2 --depend=" + str(job_id_make_fasta) + " {}/SH/full_target_mapping.sh".format(outdir))
    job_output=subprocess.getoutput("sbatch -c 2 --depend=" + str(job_id_make_fasta) + " {}/SH/full_target_mapping.sh".format(outdir))
    write_file = open("{}/SH/job_id_run_dagcon_consensus.sh".format(outdir), "w")
    write_file.write(job_output.split()[3])
    write_file.close()
