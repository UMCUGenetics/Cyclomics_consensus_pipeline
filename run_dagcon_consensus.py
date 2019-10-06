#! /usr/bin/env python
import sys, os, re
import glob
import commands
import time
import gzip
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest="wkdir", metavar="[PATH]", help="full path to FASTQ folder [default = ./")
        group.add_option("-o", dest="outdir", metavar="[PATH]", help="full path to output folder [default = ./")
        group.add_option("-c", default=10, dest="coverage", metavar="[INT]", help="minimum coverage required for assembly [default = 10]")
        group.add_option("-m", default="m.elferink@umcutrecht.nl", dest="mail", metavar="[STRING]", help="email used for job submitting [default = m.elferink@umcutrecht.nl]")
        group.add_option("-t", default="4:00:00", dest="timeslot", metavar="[TIME]", help="time slot for jobs [default = 4:00:00]")
        group.add_option("-n", default="1000", dest="number", metavar="[INT]", help="number of jobs within a scatterjob [default = 1000]")
        group.add_option("-p", dest="prefix", metavar="[STRING]", help="prefix of BAM names [default = off (FASTQ input folder name)]")
        group.add_option("--mem", default=32, dest="max_mem", metavar="[INT]", help="memory used for jobs [default = 32]")
        group.add_option("--threads", default=1, dest="threads", metavar="[INT]", help="number threads used for jobs [default = 1]")
        group.add_option("--cl", default=35, dest="cons_len", metavar="INT", help="minimum length (bp) for consensus calling [default 35]")      
        group.add_option("--project", default="compgen", dest="project", metavar="STRING", help="SGE project for submitting jobs [default compgen]")

        group.add_option("-b", default="/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa", dest="bwa", metavar="[PATH]", help="full path to bwa binary [default = /hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa ]")
        group.add_option("--sa", default="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba", dest="sambamba", metavar="[PATH]", help="full path to sambamba binary [default = /hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba]")
        group.add_option("--la", default="/hpc/cog_bioinf/ridder/tools/last-921/", dest="lastal", metavar="[PATH]", help="full path to lastal binary [default = /hpc/cog_bioinf/ridder/tools/last-921/]")

        group.add_option("-e", default="/hpc/cog_bioinf/ridder/tools/bam2m5/env_3.6/bin/activate", dest="env", metavar="[ENV]", help="full path to python enviroment [default = ihpc/cog_bioinf/ridder/tools/bam2m5_new/env_3.6/bin/activate ]")
        group.add_option("--b5", default="/hpc/cog_bioinf/ridder/tools/bam2m5/bam2m5.py", dest="bam2m5", metavar="[PATH]", help="full path to bam2m5 binary [default = /hpc/cog_bioinf/ridder/tools/bam2m5_new/bam2m5.py]")
        group.add_option("--pbdagcon", default="/hpc/cog_bioinf/ridder/tools/pbdagcon/src/cpp/pbdagcon", dest="pbdagcon", metavar="[PATH]", help="full path to pbgadcon binary [default = /hpc/cog_bioinf/ridder/tools/pbdagcon/src/cpp/pbdagcon ]")
        group.add_option("--rf", default="/hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta", dest="refgenome_full", metavar="[PATH]", help="full path to complete reference genome [default = /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome//version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta]")
        group.add_option("--rt", default="/hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/BRAF_TP53_BB_pjet.fasta", dest="refgenome_target", metavar="[PATH]", help="full path to targeted reference genome [default = /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/BRAF_TP53_BB_pjet.fasta]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

if not opt.wkdir:
    sys.exit("provide input folder")

if opt.outdir:
    outdir= opt.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
else:
    sys.exit("provide output folder")

if opt.prefix:
    runid=str(opt.prefix)
else:
    runid=wkdir.split("/")[-2]

if "fasta" not in opt.refgenome_target:
    sys.exit("please provide fasta file for target reference genome")
if "fasta" not in opt.refgenome_full:
    sys.exit("please provide fasta file for full reference genome")

wkdir=opt.wkdir
coverage= opt.coverage
number=int(opt.number)
project=opt.project
trim=0
timeslot=str(opt.timeslot)
max_mem=int(opt.max_mem)
threads=int(opt.threads)
cons_len=int(opt.cons_len)
lastal_src=str(opt.lastal)+"src/lastal"
lastparam=str(opt.lastal)+"last_params"
lastsplit=str(opt.lastal)+"src/last-split"
mafconvert=str(opt.lastal)+"scripts/maf-convert"

# Log GIT version + commit of repository
if str(sys.argv[0]) == "python":
    repo="/".join(sys.argv[1].split("/")[0:-1])
else:
    repo="/".join(sys.argv[0].split("/")[0:-1])

os.system("git --git-dir="+str(repo)+"/.git describe --tags >"+str(outdir)+"/GIT_run_dagcon_consensus.log")
os.system("git --git-dir="+str(repo)+"/.git log >>"+str(outdir)+"/GIT_run_dagcon_consensus.log")

## print log of used parameters ## 
#runid=wkdir.split("/")[-2]
write_file=open(outdir+"/"+str(runid)+".log","w")
(options, args) = parser.parse_args()
for item in vars(options):
    write_file.write(str(item)+"\t"+str(vars(options)[item])+"\n") 
write_file.close()

if wkdir.endswith('/'): # chop last "/" if present
    wkdir=wkdir[0:-1]

if os.path.isdir(outdir+"/bam") or os.path.isdir(outdir+"/m5") or os.path.isdir(outdir+"/consensus"):
    sys.exit("please remove previous analysis first")

refgenome_target_db=opt.refgenome_target[0:-6]
refgenome_target=opt.refgenome_target
list=[]

y=1
runner=0
files= commands.getoutput("find "+wkdir+ " -iname \"*fastq\"").split()
if len(files)==0: ## likely that fastq are compressed
    files= commands.getoutput("find "+wkdir+ " -iname \"*fastq.gz\"").split()
    gz="on"    
else:
    gz="off"

fastq_dic={}
line_dic={}
for f in files:
    if gz=="off":
        lines=open(f, "r").readlines()
    else:
        lines=gzip.open(f, "rb").readlines()
    x=0
    while x< len(lines):
        if "@" in lines[x]:	# bypass nasty bug in which fastq's are corrupted. Corrupted will be skipped, not fixed
            splitline=lines[x].split()
            read=str(splitline[0].split("@")[1])+"_read_"+splitline[3].split("=")[1]+"_ch_"+splitline[4].split("=")[1]
            if str(f) not in line_dic:
                line_dic[str(f)]=[{read:[x+1,x+4]}]
            else:
                line_dic[str(f)]+=[{read:[x+1,x+4]}]
            x+=4
        else:
            print "Error ",f,x
            y-=1
            x+=1
        y+=1
    runner+=1
total =0

for item in line_dic:
    if gz=="off":
        if "pass" in item or "fail" in item:  # handles run with or without pass/fail reads
            folder=item.split("/")[-2] +"_"+item.split("/")[-1][0:-6]
        else:
            folder=item.split("/")[-1][0:-6]
    elif gz=="on":
        if "pass" in item or "fail" in item: # handles run with or without pass/fail reads
            folder=item.split("/")[-2] +"_"+item.split("/")[-1][0:-9]
        else:
            folder=item.split("/")[-1][0:-9]
 
    x=1
    c=0
    write_file=open(outdir+"/x"+folder+"_job_"+str(x)+".sh","w")
    write_file.write(". "+opt.env+"\n")
    write_file.write("echo \"Start poststats    \" `date` \"    \" `uname -n`\n\n")
    os.system("mkdir -p "+ outdir+"/bam/"+ folder+"_"+str(x))
    os.system("mkdir -p "+ outdir+"/m5/"+ folder+"_"+str(x))
    os.system("mkdir -p "+ outdir+"/consensus/"+ folder+"_"+str(x))

    for f in line_dic[item]:
        total+=1
        position=f.values()[0]
        fastq=f.keys()[0]

        ### Make BAM file
        if gz=="off":
            write_file.write("cat "+str(item)+"| sed -n \'"+str(position[0])+","+ str(position[1])+"\'p | "+ lastal_src+" -Q 1 -p "+lastparam+" "+refgenome_target_db + " -P 1 /dev/stdin |"+ lastsplit+ "|" + mafconvert+ " -f "+ refgenome_target_db +".dict sam -r \"ID:"+ fastq +" PL:nanopore SM:"+ fastq +"\" /dev/stdin | "+ opt.sambamba+ " view -S -f bam /dev/stdin | "+ opt.sambamba + " sort -t 1 /dev/stdin -o "+ outdir+"/bam/"+ folder+"_"+str(x) +"/"+fastq+".sorted.bam\n")
        elif gz=="on":
           write_file.write("zcat "+str(item)+"| sed -n \'"+str(position[0])+","+ str(position[1])+"\'p | "+ lastal_src+" -Q 1 -p "+lastparam+" "+refgenome_target_db + " -P 1 /dev/stdin |"+ lastsplit+ "|" + mafconvert+ " -f "+ refgenome_target_db +".dict sam -r \"ID:"+ fastq +" PL:nanopore SM:"+ fastq +"\" /dev/stdin | "+ opt.sambamba+ " view -S -f bam /dev/stdin | "+ opt.sambamba + " sort -t 1 /dev/stdin -o "+ outdir+"/bam/"+ folder+"_"+str(x) +"/"+fastq+".sorted.bam\n") 

        ## Add Bam to TAR ball
        write_file.write("cd "+outdir+"/bam/"+ folder+"_"+str(x) +"/\n") 
        write_file.write("tar  --remove-files -f "+folder+"_"+str(x)+".tar"+ " -r "+fastq+".sorted.bam*\n")

        # Extract BAM file from TAR ball and stdout to command to make m5
        write_file.write("tar -axf "+outdir+"/bam/"+ folder+"_"+str(x) +"/"+folder+"_"+str(x)+".tar "+ fastq+".sorted.bam -O | python "+opt.bam2m5+ " - "+refgenome_target+ " "+ outdir+"/m5/"+ folder+"_"+str(x) +"/"+ fastq+".sorted.m5\n")

        ## Add m5 to TAR ball
        write_file.write("cd "+outdir+"/m5/"+ folder+"_"+str(x) +"/\n")
        write_file.write("tar  --remove-files -f "+folder+"_"+str(x)+".tar"+ " -r "+fastq+".sorted.m5*\n")

        # Extract M5 file from TAR ball and stdout to command
        write_file.write("tar -axf "+outdir+"/m5/"+ folder+"_"+str(x) +"/"+folder+"_"+str(x)+".tar "+ fastq+".sorted.m5 -O |" + opt.pbdagcon+ " - "+" -m "+str(cons_len)+" -c "+str(coverage)+" -t "+str(trim)+" -j "+str(threads)+" > " + outdir+"/consensus/"+ folder+"_"+str(x) +"/"+ fastq+".consensus\n")
        write_file.write("sed -i \'s/>/>"+f.keys()[0]+"_/g\' "+ outdir+"/consensus/"+ folder+"_"+str(x) +"/"+fastq+".consensus\n")

        ## Add consensus to TAR ball
        write_file.write("cd "+outdir+"/consensus/"+ folder+"_"+str(x) +"/\n")
        write_file.write("tar  --remove-files -f "+folder+"_"+str(x)+".tar"+ " -r "+fastq+".consensus*\n")

        if c==number-1:
            list+=[outdir+"/x"+folder+"_job_"+str(x)+".sh"]
            write_file.close()
            x+=1
            if total == len(line_dic[item]):
                pass
            else:
                write_file=open(outdir+"/x"+folder+"_job_"+str(x)+".sh","w")
                write_file.write(". "+opt.env+"\n")
                write_file.write("echo \"Start poststats    \" `date` \"    \" `uname -n`\n")
                os.system("mkdir "+ outdir+"/bam/"+ folder+"_"+str(x))
                os.system("mkdir "+ outdir+"/m5/"+ folder+"_"+str(x))
                os.system("mkdir "+ outdir+"/consensus/"+ folder+"_"+str(x))

            c=0
        else:
            c+=1

    list+=[outdir+"/x"+folder+"_job_"+str(x)+".sh"]
    write_file.close()

os.mkdir(outdir+"/SH")
job_id=[]
for job in list:
    action= "qsub -q all.q -P "+ str(project)+ " -l h_rt=" +str(timeslot)+ " -l h_vmem="+str(max_mem)+"G -R y -cwd -pe threaded "+str(threads)+" "+str(job) + " -o "  +str(outdir+"/SH/")+ " -e " + str(outdir+"/SH/" + " -m a -M "+str(opt.mail))
    job_id+= [commands.getoutput(action).split()[2]]

#runid=wkdir.split("/")[-1]
write_file=open(outdir+"/full_target_mapping.sh","w")
write_file.write("mkdir "+ str(outdir)+"/tmp\n")

# extract all file from all consensus TAR balls and put in FASTA file
write_file.write("find "+outdir +"/consensus/ -type f -iname \"*tar\" -exec tar -O -xf {} \; >> "+outdir+"/"+runid+"_full_consensus.fasta\n")
write_file.write(opt.bwa + " mem -t 2 -c 100 -M -R \"@RG\\tID:"+runid+"\\tSM:"+runid+"\\tPL:NANOPORE\\tLB:"+runid+"\" "+opt.refgenome_full+" "+outdir+"/"+runid+"_full_consensus.fasta | " + opt.sambamba+ " view -S -f bam /dev/stdin |"+  opt.sambamba+ " sort -t 2 --tmpdir="+str(outdir)+"/tmp /dev/stdin -o "+outdir+"/"+runid+"_full_consensus.sorted.bam\n")
write_file.write("gzip "+outdir+"/"+runid+"_full_consensus.fasta\n") 
write_file.close()
action=("qsub -cwd -q all.q -P "+str(project)+ " -l h_rt=36:00:00 -l h_vmem=40G "+" -hold_jid "+str(",".join(job_id))+" "+str(outdir)+"/full_target_mapping.sh" + " -o "  +str(outdir)+"/SH/"+ " -e " + str(outdir)+"/SH/" + " -m ae -M "+str(opt.mail))
hold_id=[commands.getoutput(action).split()[2]]

write_file=open(outdir+"/cleanup.sh","w")
write_file.write("mv "+ str(outdir)+"/*sh* "+str(outdir)+"/SH\n")
write_file.write("cd "+str(outdir)+"/SH\n") 
write_file.write("zip -m SH.zip *\n")
write_file.close()

action=("qsub -cwd -q all.q -P "+str(project)+ " -l h_rt=0:05:00 -l h_vmem=1G "+" -hold_jid "+str(",".join(hold_id))+" "+str(outdir)+"/cleanup.sh" + " -o "  +str(outdir)+"/SH/"+ " -e " + str(outdir)+"/SH/" + " -m baes -M "+str(opt.mail))
commands.getoutput(action)

print "All jobs sumbitted" 



