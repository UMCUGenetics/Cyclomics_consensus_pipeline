#! /usr/bin/env python
import sys, os,re, commands
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest="wkdir", metavar="[PATH]", help="full path to  m5 folder [default = ./]")
        group.add_option("-o", dest="outdir", metavar="[PATH]", help="full path to output folder [default = ./")
        group.add_option("-p", default="/hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon", dest="pbdagcon", metavar="[PATH]", help="full path to pbgadcon binary [default = /hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon]")
        group.add_option("--param", default=" -c 1 -t 0 -j 2 ", dest="param", metavar="[STRING]", help="pbdagcon paramteres [default =  -c 2 -t 0 -j 2 ]")
        group.add_option("--sa", default="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba", dest="sambamba", metavar="[PATH]", help="full path to sambamba binary [default = /hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba]")
        group.add_option("-b", default="/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa", dest="bwa", metavar="[PATH]", help="full path to bwa binary [default = /hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa ]")
        group.add_option("--rf", default="/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta", dest="refgenome_full", metavar="[PATH]", help="full path to complete reference genome [default = /hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta]")
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
param=" -m "+str(opt.cons_len)+str(opt.param)

if opt.wkdir:
    cwd=opt.wkdir
else:
    cwd = os.getcwd()

if opt.outdir:
    outfolder=opt.outdir
else:
    outfolder = os.getcwd()

infolder = opt.wkdir
sambamba=opt.sambamba
bwa=opt.bwa
refgenome_full=opt.refgenome_full
os.system("mkdir "+str(outfolder)+"/tmp")
os.system("mkdir "+str(outfolder)+"/SH")
os.system("mkdir "+str(outfolder)+"/bin_consensus")
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


for folder in os.listdir(infolder):
    tarfile = str(infolder)+"/"+str(folder)+"/"+str(folder).split("/")[-1]+".tar"
    files=commands.getoutput("tar -tf "+str(tarfile)).split()
    for f in files:
        in_file=commands.getoutput("tar -axf "+str(tarfile)+ " "+ str(f)+ " -O ").split()
        insert_count=0
        for line in in_file:
            if str(insert) in line: 
                ## count specific reads for variable 'insert'. e.g. if TP53, backbone counts will not be taken into account
                ## pitfall for multi-exon experiments: e.g. all TP53 insert will be counted as insert. Could be improved if mapping location if also included.
                insert_count+=1
        ## notice that in pbdagcon both insert and backbone are parsed into fasta. Only one 1 is count-selection performed
        if insert_count == 0:
            pass
        elif insert_count >= 40:
            os.system("tar -axf "+str(tarfile)+ " "+ str(f)+ " -O | " + str(pbdagcon)+" - "+ param + " | sed 's/>/>"+str(f[0:-10])+"_"+"/g' 1>> "+str(outfolder)+"//bin_consensus/consensus_40+.fasta")
        else:
            for x in xrange(1, 40):
                if insert_count == x:
                    os.system("tar -axf "+str(tarfile)+ " "+ str(f)+ " -O | " + str(pbdagcon)+" - "+ param + " | sed 's/>/>"+str(f[0:-10])+"_"+"/g' 1>> "+str(outfolder)+"//bin_consensus/consensus_"+str(x)+".fasta")

job_id=[]
def write_new_file(x,job_id):
    runid="consensus_"+str(x)
    write_file=open(str(outfolder)+"/SH/"+str(runid)+".sh","w")
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
    action="qsub -q all.q -P "+str(project)+" -l h_rt="+str(opt.timeslot)+" -l h_vmem="+str(opt.max_mem)+"G -cwd -pe threaded "+str(opt.threads)+" "+str(outfolder)+"/SH/"+str(runid+".sh") + " -o "+str(outfolder)+"/SH/"+str(runid+".output")+ " -e ./SH/" + str(runid+".error") +" -m a -M "+ str(mail)
    job_id+= [commands.getoutput(action).split()[2]]
    return job_id

for x in xrange(1, 40):
    job_id=write_new_file(x,job_id)
job_id=write_new_file("40+",job_id)

# cleanup
write_file=open(str(outfolder)+"/SH/cleanup.sh","w")
write_file.write("mv "+str(outfolder)+"/*sh* SH\n")
write_file.close()
action=("qsub -cwd -q all.q -P "+str(project)+ " -l h_rt=00:05:00 -l h_vmem=4G "+" -hold_jid "+str(",".join(job_id))+" "+str(outfolder)+"/SH/cleanup.sh" + " -o "  +str(outfolder)+"/SH/"+ " -e " + str(outfolder)+"/SH/" + " -m ae -M "+str(opt.mail))
os.system(action)
