#! /usr/bin/env python
import sys, os,re, commands
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("--im", dest="m5dir", metavar="[PATH]", help="full path to  m5 folder [default = ./]")
        group.add_option("--ib", dest="bamdir", metavar="[PATH]", help="full path to  bam folder [default = ./]")
        group.add_option("-o", dest="outdir", metavar="[PATH]", help="full path to output folder [default = ./")
        group.add_option("-p", default="/hpc/cog_bioinf/ridder/tools/pbdagcon/src/cpp/pbdagcon", dest="pbdagcon", metavar="[PATH]", help="full path to pbgadcon binary [default = /hpc/cog_bioinf/ridder/tools/pbdagcon/src/cpp/pbdagcon]")
        #group.add_option("--param", default=" -c 1 -t 0 -j 2 ", dest="param", metavar="[STRING]", help="pbdagcon paramteres [default =  -c 2 -t 0 -j 2 ]")
        group.add_option("--sa", default="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba", dest="sambamba", metavar="[PATH]", help="full path to sambamba binary [default = /hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba]")
        group.add_option("-b", default="/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa", dest="bwa", metavar="[PATH]", help="full path to bwa binary [default = /hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa ]")
        group.add_option("--rf", default="/hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta", dest="refgenome_full", metavar="[PATH]", help="full path to complete reference genome [default = /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta]")
        group.add_option("--project", default="compgen", dest="project", metavar="STRING", help="SGE project for submitting jobs [default compgen]")
        group.add_option("-m", default="m.elferink@umcutrecht.nl", dest="mail", metavar="[STRING]", help="email used for job submitting [default = m.elferink@umcutrecht.nl]")
        group.add_option("-a", default="TP53", dest="insert", metavar="[STRING]", help="name of insert chromosome [default = TP53]")
        group.add_option("-t", default="4:00:00", dest="timeslot", metavar="[TIME]", help="time slot for jobs [default = 4:00:00]")
	group.add_option("--mem", default=32, dest="max_mem", metavar="[INT]", help="memory used for jobs [default = 32]")
	group.add_option("--threads", default=4, dest="threads", metavar="[INT]", help="number threads used for jobs [default = 4]")
        group.add_option("--cl", default=35, dest="cons_len", metavar="INT", help="minimum length (bp) for consensus calling [default 35]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()


pbdagcon=opt.pbdagcon
#param=" -m "+str(opt.cons_len)+str(opt.param)


param1= " -m "+str(opt.cons_len)+" -c "
param2= " -t 0 -j 2 "

bamfolder=opt.bamdir

m5folder=opt.m5dir
bamfolder=opt.bamdir

#if opt.wkdir:
#    cwd=opt.wkdir
#else:
#    cwd = os.getcwd()

if opt.outdir:
    outfolder=opt.outdir
else:
    outfolder = os.getcwd()

sambamba=opt.sambamba
bwa=opt.bwa
refgenome_full=opt.refgenome_full
os.system("mkdir "+str(outfolder)+"/tmp")
os.system("mkdir "+str(outfolder)+"/SH")
os.system("mkdir "+str(outfolder)+"/bin_consensus")
os.system("mkdir "+str(outfolder)+"/bin_consensus_folder")

mail=opt.mail
project=opt.project
insert=opt.insert


# Log GIT version + commit of repository
if str(sys.argv[0]) == "python":
    repo="/".join(sys.argv[1].split("/")[0:-1])
else:
    repo="/".join(sys.argv[0].split("/")[0:-1])

os.system("git --git-dir="+str(repo)+"/.git describe --tags >"+str(outfolder)+"/GIT_bin_on_repeat_count.log")
os.system("git --git-dir="+str(repo)+"/.git log >>"+str(outfolder)+"/GIT_bin_on_repeat_count.log")

job_id=[]
for folder in os.listdir(bamfolder):
    if os.path.isfile(str(bamfolder)+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"):
        bamtarfile = str(bamfolder)+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"
        m5tarfile = str(m5folder)+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"
        files=commands.getoutput("tar -tf "+str(bamtarfile)).split()
        write_file=open(str(outfolder)+"/SH/"+str(folder)+".sh","w")
        for f in files:
            insert_count=0
            if "bai" not in str(f):
                os.system("tar -axf "+str(bamtarfile)+ " "+ str(f)+"*")
                insert_count=commands.getoutput(str(sambamba)+ " depth base -L "+str(insert)+ " " + str(f) + "| cut -f3| sort -nk3 | head -n1").split("\n")
                if len(insert_count)> 1 and "failed" not in str(insert_count):
                    insert_count=insert_count[1]
                else:
                    insert_count=0  
                os.system("rm "+ str(f)+"*")
            
                if str(insert_count) == "0" or str(insert_count) == "COV":
                    pass
                elif int(insert_count) >= 40:
                    threshold = 38
                    write_file.write("tar -axf "+str(m5tarfile)+ " "+ str(f[0:-11])+ "* -O | " + str(pbdagcon)+" - "+ str(param1)+ str(threshold)+ str(param2) + " | sed 's/>/>"+str(f[0:-10])+"_"+"/g' 1>> "+str(outfolder)+"/bin_consensus_folder/"+str(folder)+"_consensus_40+.fasta\n")
                else:
                    for x in xrange(1, 40):
                        if x > 3:
                            threshold = x - 2
                        else:
                            threshold = 1
                        if int(insert_count) == x:
                            write_file.write("tar -axf "+str(m5tarfile)+ " "+ str(f[0:-11])+ "* -O | " + str(pbdagcon)+" - "+ str(param1)+ str(threshold)+ str(param2) + " | sed 's/>/>"+str(f[0:-10])+"_"+"/g' 1>> "+str(outfolder)+"/bin_consensus_folder/"+str(folder)+"_consensus_"+str(x)+".fasta\n")

    write_file.close()
    action="qsub -q all.q -P "+str(project)+" -l h_rt="+str(opt.timeslot)+" -l h_vmem=10G -cwd -pe threaded 1 -o "+str(outfolder)+"/SH/"+str(folder)+".output -e "+ str(outfolder)+"/SH/"+str(folder)+".error -m a -M "+ str(mail) + " "+ str(outfolder)+"/SH/"+str(folder)+".sh"
    job_id+= [commands.getoutput(action).split()[2]]



write_file=open(str(outfolder)+"/SH/merge_fasta.sh","w")
for x in xrange(1, 40):
    write_file.write("find "+str(outfolder)+"/bin_consensus_folder/ -iname \"*_consensus_"+str(x)+".fasta\" -exec cat {} \; >> "+str(outfolder)+"/bin_consensus/consensus_"+str(x)+".fasta\n")
write_file.write("find "+str(outfolder)+"/bin_consensus_folder/ -iname \"*_consensus_40+.fasta\" -exec cat {} \; >> "+str(outfolder)+"/bin_consensus/consensus_40+.fasta\n")
write_file.close()
action= "qsub -q all.q -P "+str(project)+" -l h_rt=2:0:0 -l h_vmem=20G -cwd -pe threaded 1 -o "+ str(outfolder)+"/SH/merge_fasta.output -e "+ str(outfolder)+"/SH/merge_fasta.error -m a -M "+ str(mail) + " -hold_jid "+str(",".join(job_id))+" "+str(outfolder)+"/SH/merge_fasta.sh" 
job_id_merge=commands.getoutput(action).split()[2]

job_id=[]
def write_new_file(x,job_id):
    runid="consensus_"+str(x)
    write_file=open(str(outfolder)+"/SH/"+str(runid)+"_mapping.sh","w")
    write_file.write(bwa + " mem -t "+str(opt.threads)+" -c 100 -M -R \"@RG\\tID:"+runid+"\\tSM:"+runid+"\\tPL:NANOPORE\\tLB:"+runid+"\" "+refgenome_full+" "+str(outfolder)+"//bin_consensus/"+runid+".fasta > "+str(outfolder)+"/"+runid+"_full_consensus.sam\n")
    write_file.write(sambamba+ " view -S -f bam "+str(outfolder)+"/" +runid+"_full_consensus.sam > "+str(outfolder)+"/"+runid+"_full_consensus.bam\n")
    write_file.write(sambamba+ " sort -t "+str(opt.threads)+" --tmpdir=./tmp"+" "+str(outfolder)+"/"+runid+"_full_consensus.bam -o "+str(outfolder)+"/"+runid+"_full_consensus.sorted.bam\n")
    write_file.write("sleep 2\n")
    write_file.write("rm "+str(outfolder)+"/"+runid+"_full_consensus.sam\n")
    write_file.write("sleep 2\n")
    write_file.write("rm "+str(outfolder)+"/"+runid+"_full_consensus.bam\n")
    write_file.write("sleep 2\n")
    write_file.write("mv "+str(outfolder)+"/"+runid+"_full_consensus.sorted.bam* "+str(outfolder)+"/bin_consensus/\n")
    write_file.write("sleep 2\n")
    write_file.close()
    action="qsub -q all.q -P "+str(project)+" -l h_rt="+str(opt.timeslot)+" -l h_vmem="+str(opt.max_mem)+"G -cwd -pe threaded "+str(opt.threads)+ " -o "+str(outfolder)+"/SH/"+str(runid)+".output"+ " -e "+str(outfolder)+"/SH/" + str(runid)+".error" +" -m a -M "+ str(mail)+ " -hold_jid "+str(job_id_merge) + " "+str(outfolder)+"/SH/"+str(runid)+"_mapping.sh"
    #print "jrrp", action
    job_id+= [commands.getoutput(action).split()[2]]
    return job_id

for x in xrange(1, 40):
    #if os.path.isfile(str(outfolder)+"/bin_consensus/consensus_"+str(x)+".fasta"):
    job_id=write_new_file(x,job_id)

#if os.path.isfile(str(outfolder)+"/bin_consensus/consensus_40+.fasta"):
#    job_id=write_new_file("40+",job_id)

job_id=write_new_file("40+",job_id)


# cleanup
write_file=open(str(outfolder)+"/SH/cleanup.sh","w")
write_file.write("mv "+str(outfolder)+"/*sh* SH\n")
write_file.close()
action=("qsub -cwd -q all.q -P "+str(project)+ " -l h_rt=00:05:00 -l h_vmem=4G "+" -hold_jid "+str(",".join(job_id))+" -o "  +str(outfolder)+"/SH/"+ " -e " + str(outfolder)+"/SH/" + " -m baes -M "+str(opt.mail)) + " "+str(outfolder)+"/SH/cleanup.sh"
os.system(action)
