#! /usr/bin/env python
import sys, os, re
import glob
#import commands
import subprocess
import time
import gzip
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
    parser = OptionParser();
    group = OptionGroup(parser, "Main options")
    group.add_option("-i", dest = "wkdir", metavar = "[PATH]", help = "full path to FASTQ folder [default = ./")
    group.add_option("-o", dest = "outdir", metavar = "[PATH]", help = "full path to output folder [default = ./")
    group.add_option("-b", dest = "blacklist", metavar = "[PATH]", help = "full path to blacklist file [default = off]")
    group.add_option("-c", default = 10, dest = "coverage", metavar = "[INT]", help = "minimum coverage required for assembly [default = 10]")
    group.add_option("-m", default = "m.elferink@umcutrecht.nl", dest = "mail", metavar = "[STRING]", help = "email used for job submitting [default = m.elferink@umcutrecht.nl]")
    #group.add_option("-t", default = "4:00:00", dest = "timeslot", metavar = "[TIME]", help = "time slot for jobs [default = 4:00:00]")
    group.add_option("-n", default = "1000", dest = "number", metavar="[INT]", help = "number of jobs within a scatterjob [default = 1000]")
    group.add_option("-p", dest = "prefix", metavar = "[STRING]", help = "prefix of BAM names [default = off (FASTQ input folder name)]")
    #group.add_option("--mem", default = 32, dest = "max_mem", metavar = "[INT]", help = "memory used for jobs [default = 32]")
    group.add_option("--threads", default = 1, dest = "threads", metavar = "[INT]", help = "number threads used for jobs [default = 1]")
    group.add_option("--cl", default = 35, dest = "cons_len", metavar = "INT", help = "minimum length (bp) for consensus calling [default 35]")      
    #group.add_option("--project", default = "compgen", dest = "project", metavar = "STRING", help = "SGE project for submitting jobs [default compgen]")
    group.add_option("--bwa", default = "/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa", dest = "bwa", metavar = "[PATH]", help = "full path to bwa binary [default = /hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa ]")
    group.add_option("--sa", default = "/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba", dest = "sambamba", metavar = "[PATH]", help = "full path to sambamba binary [default = /hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba]")
    group.add_option("--la", default = "/hpc/compgen/tools/last-921/", dest = "lastal", metavar = "[PATH]", help = "full path to lastal binary [default = /hpc/compgen/tools/last-921/]")
    #group.add_option("-e", default = "/hpc/compgen/tools/bam2m5/env_3.6/bin/activate", dest = "env", metavar = "[ENV]", help = "full path to python enviroment [default = ihpc/compgen/tools/bam2m5_new/env_3.6/bin/activate ]")
    group.add_option("--b5", default = "/hpc/compgen/tools/bam2m5/bam2m5.py", dest = "bam2m5", metavar = "[PATH]", help = "full path to bam2m5 binary [default = /hpc/compgen/tools/bam2m5_new/bam2m5.py]")
    group.add_option("--pbdagcon", default = "/hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon", dest = "pbdagcon", metavar = "[PATH]", help = "full path to pbgadcon binary [default = /hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon ]")
    group.add_option("--rf", default = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta", dest = "refgenome_full", metavar = "[PATH]", help = "full path to complete reference genome [default = /hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta]")
    group.add_option("--rt", default = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta", dest = "refgenome_target", metavar = "[PATH]", help = "full path to targeted reference genome [default = /hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta]")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

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
                    print("Warning, read found multiple times in blacklist {}".format(item.split()[1].split("_")[1]))

    coverage = opt.coverage
    number = int(opt.number)
    #project = opt.project
    trim = 0
    #timeslot = str(opt.timeslot)
    #max_mem = int(opt.max_mem)
    threads = int(opt.threads)
    cons_len = int(opt.cons_len)
    lastal_src = str(opt.lastal) + "src/lastal"
    lastparam = str(opt.lastal) + "last_params"
    lastsplit = str(opt.lastal) + "src/last-split"
    mafconvert = str(opt.lastal) + "scripts/maf-convert"

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
        files =subprocess.getoutput("find {} -iname \"*fastq.gz\"".format(wkdir)).split()
        gz = "on"    
    else:
        gz = "off"

    fastq_dic = {}
    line_dic = {}
    #job_list = []
    linecount = 1
    runner = 0
    for f in files:
        if gz == "off":
            lines = open(f, "r").readlines()
        else:
            lines = gzip.open(f, "rb").readlines()
        fastqline = 0
        while fastqline < len(lines):
            if "@" in lines[fastqline].decode("utf-8"):	# bypass nasty bug in which fastq's are corrupted. Corrupted will be skipped, not fixed
                splitline = lines[fastqline].decode("utf-8").split()
                read = "{0}_read_{1}_ch_{2}".format(splitline[0].split("@")[1], splitline[3].split("=")[1], splitline[4].split("=")[1])
                if f not in line_dic:
                    line_dic[f] = [{read:[fastqline + 1, fastqline + 4]}]
                else:
                    line_dic[f] += [{read:[fastqline + 1, fastqline + 4]}]
                fastqline += 4
            else:
                print("Error ", f, fastqline)
                linecount -= 1
                fastqline += 1
            linecount += 1
        runner += 1

    total = 0
    for item in line_dic:
        if gz == "off":
            if "pass" in item or "fail" in item:  # handles run with or without pass/fail reads
                folder = "{0}_{1}".format(item.split("/")[-2], item.split("/")[-1][0:-6])
            else:
                folder = item.split("/")[-1][0:-6]
        elif gz == "on":
            if "pass" in item or "fail" in item: # handles run with or without pass/fail reads
                folder = "{0}_{1}".format(item.split("/")[-2], item.split("/")[-1][0:-9])
            else:
                folder = item.split("/")[-1][0:-9]
 
        subnumber = 1
        count = 0
        #write_file = open("{0}/x{1}_job_{2}.sh".format(outdir, folder, subnumber), "w")
        #write_file.write(". {}\n".format(opt.env))
        #write_file.write("echo \"Start poststats    \" `date` \"    \" `uname -n`\n\n")
        os.system("mkdir -p {0}/bam/{1}_{2}".format(outdir, folder, subnumber))
        os.system("mkdir -p {0}/m5/{1}_{2}".format(outdir, folder, subnumber))
        os.system("mkdir -p {0}/consensus/{1}_{2}".format(outdir, folder, subnumber))
        for f in line_dic[item]:
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
                action = ("{cat} {item} | sed -n \'{pos0},{pos1}\'p | {lastall} -Q 1 -p {lastparam} {refgenome_target_db} -P 1 /dev/stdin | {lastsplit} | {mafconvert} -f {refgenome_target_db}.dict sam -r \"ID:{fastq} PL:nanopore SM:{fastq}\" /dev/stdin | {sambamba} view -S -f bam /dev/stdin | {sambamba} sort -t 1 /dev/stdin -o {outdir}/bam/{folder}_{subnumber}/{fastq}.sorted.bam ".format(
                        cat = cat,
                        item = item,
                        pos0 = position[0],
                        pos1 = position[1],
                        lastall = lastal_src,
                        lastparam = lastparam,
                        refgenome_target_db = refgenome_target_db,
                        lastsplit = lastsplit,
                        mafconvert = mafconvert,
                        fastq = fastq,
                        sambamba = opt.sambamba,
                        outdir = outdir,
                        folder = folder,
                        subnumber = subnumber,
                        ))
                os.system(action)
                #print(action)

                """Add Bam to TAR ball"""
                
                #action = ("cd {outdir}/bam/{folder}_{subnumber}\n".format(
                #    outdir = outdir,
                #    folder = folder,
                #    subnumber = subnumber
                #    ))

                os.chdir("{outdir}/bam/{folder}_{subnumber}".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber
                    ))

                #os.system(action)
                #print(action)

                action =("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.sorted.bam*".format(
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq
                    ))
                os.system(action)
                #print(action)                 

                """Extract BAM file from TAR ball and stdout to command to make m5"""
                action = ("tar -axf {outdir}/bam/{folder}_{subnumber}/{folder}_{subnumber}.tar {fastq}.sorted.bam -O | python {bam2m5} - {refgenome_target} {outdir}/m5/{folder}_{subnumber}/{fastq}.sorted.m5".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq,
                    bam2m5 = opt.bam2m5,
                    refgenome_target = refgenome_target
                    ))
                os.system(action) 
                #print(action)

                """Add m5 to TAR ball"""
                #action = ("cd {outdir}/m5/{folder}_{subnumber}\n".format(
                #    outdir = outdir,
                #    folder = folder,
                #    subnumber = subnumber
                #    ))
                #os.system(action) 
                #print(action)
                os.chdir("{outdir}/m5/{folder}_{subnumber}".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber
                    ))

                action = ("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.sorted.m5*".format(
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq
                    ))
                os.system(action) 
                #print(action) 
 
                """Extract M5 file from TAR ball and stdout to command"""
                action = ("echo {fastq} >> {outdir}/consensus/Full_insert_count_{folder}_{subnumber}.txt".format(
                    fastq = fastq,
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber
                    ))
                os.system(action)
                #print(action)

                action = ("tar -axf {outdir}/m5/{folder}_{subnumber}/{folder}_{subnumber}.tar {fastq}.sorted.m5 -O |{pbdagcon} - -m {cons_len} -c {coverage} -t {trim} -j {threads} 1> {outdir}/consensus/{folder}_{subnumber}/{fastq}.consensus 2>> {outdir}/consensus/Full_insert_count_{folder}_{subnumber}.txt".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq,
                    pbdagcon = opt.pbdagcon,
                    cons_len = cons_len,
                    coverage = coverage,
                    trim = trim,
                    threads = threads
                    ))
                os.system(action)
                #print(action)

                action = ("sed -i \'s/>/>{fastq}_/g\' {outdir}/consensus/{folder}_{subnumber}/{fastq}.consensus".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq
                    ))
                os.system(action)
                #print(action)    

                """Add consensus to TAR ball"""
                #action = ("cd {outdir}/consensus/{folder}_{subnumber}\n".format(
                #    outdir = outdir,
                #    folder = folder,
                #    subnumber = subnumber
                #    ))
                #os.system(action)
                #print(action)
                os.chdir("{outdir}/consensus/{folder}_{subnumber}".format(
                    outdir = outdir,
                    folder = folder,
                    subnumber = subnumber
                    ))


                action = ("tar --remove-files -f {folder}_{subnumber}.tar -r {fastq}.consensus*".format(
                    folder = folder,
                    subnumber = subnumber,
                    fastq = fastq
                    ))
                os.system(action)  
                #print(action)
          

            if count == (number - 1):
                subnumber += 1
                if total == len(line_dic[item]):
                    pass
                else:
                    os.system("mkdir -p {0}/bam/{1}_{2}".format(outdir, folder, subnumber))
                    os.system("mkdir -p {0}/m5/{1}_{2}".format(outdir, folder, subnumber))
                    os.system("mkdir -p {0}/consensus/{1}_{2}".format(outdir, folder, subnumber))
                count = 0
            else:
                count += 1

    action = ("mkdir {}/tmp\n".format(outdir))
    os.system(action)

    """Extract all file from all consensus TAR balls and put in FASTA file"""
    action = ("find {outdir}/consensus/ -type f -iname \"*tar\" -exec tar -O -xf {{}} \; >> {outdir}/{runid}_full_consensus.fasta".format(
        runid = runid,
        outdir = outdir
        ))  
    os.system(action)

    action = ("{bwa} mem -t 2 -c 100 -M -R \"@RG\\tID:{runid}\\tSM:{runid}\\tPL:NANOPORE\\tLB:{runid}\" {refgenome_full} {outdir}/{runid}_full_consensus.fasta | {sambamba} view -S -f bam /dev/stdin | {sambamba} sort -t 2 --tmpdir={outdir}/tmp /dev/stdin -o {outdir}/{runid}_full_consensus.sorted.bam".format(
        bwa = opt.bwa,
        runid = runid,
        refgenome_full = opt.refgenome_full,
        outdir = outdir,
        sambamba = opt.sambamba
        ))
    os.system(action)
    os.system("gzip {outdir}/{runid}_full_consensus.fasta\n".format(outdir = outdir, runid = runid))

    print("All jobs sumbitted")
