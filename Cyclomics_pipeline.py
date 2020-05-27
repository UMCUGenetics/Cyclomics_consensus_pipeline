#! /usr/bin/env python
import sys
import os
import subprocess

cwd = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/"
timeslot1 = "4:0:0"
timeslot2 = "8:0:0"
timeslot3 = "78:0:0"
project = "compgen"
mem = 20
cores = 2 
overlap = 100
number_reads = 100000 # must be less than total readcount


mail = "m.elferink@umcutrecht.nl"
venv = str(cwd) + "/venv/bin/activate"

default = str(cwd) + "/run_dagcon_consensus.py"
split = str(cwd) + "/split_forward_reverse_reads.py"
repeat = str(cwd) + "/bin_on_repeat_count.py"
calculate = str(cwd) + "/calculate_depth.py"
plot_dashboard = str(cwd) + "/plot_Dashboard.R"
check_numbers = str(cwd) + "/check_numbers.py"
find_read_bam = str(cwd) + "/find_read_bam.py"
structure = str(cwd) + "/make_structure.py"

rscript= "/hpc/local/CentOS7/common/lang/R/3.2.2/bin/Rscript"
full_reference = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta"
target_reference = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta"

input_folder = sys.argv[1]
output_folder = sys.argv[2]
prefix =  sys.argv[3]
insert_locus = sys.argv[4]  #e.g. TP53
backbone_locus = sys.argv[5]  #e.g. BB25
insert_targetinterval = sys.argv[6] #e.g. TP53:1-27760 > targeted reference genome!

"""Default consensus calling"""
action = "source {venv} && {default} -i {input_folder} -o {output_folder} -p {prefix} --rf={full_reference} --rt={target_reference}".format(
    venv=venv,
    default=default,
    input_folder=input_folder,
    output_folder=output_folder,
    prefix=prefix,
    full_reference=full_reference,
    target_reference=target_reference
)
os.system(action)
jobid_default = open("{output_folder}/SH/job_id_run_dagcon_consensus.sh".format(output_folder=output_folder), "r").readline().strip()
print("Submitted default consensus calling\n",action,jobid_default)


""" Forward and Reverse splitting """
""" Insert """
os.system("mkdir {output_folder}/for_rev_split_insert".format(output_folder=output_folder))
os.chdir("{output_folder}/for_rev_split_insert".format(output_folder=output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {cores} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_insert/ --rf={full_reference} -g {insert_locus} \"".format(
    timeslot=timeslot2,
    project=project,
    mem=mem,
    cores=cores,
    mail=mail,
    depend=jobid_default,
    venv=venv,
    split=split,
    output_folder=output_folder,
    full_reference=full_reference,
    insert_locus=insert_locus
)
jobid =subprocess.getoutput(action)
jobid_split_insert = jobid.split()[3]
print("Submitted forward reverse split insert\n",action,jobid_split_insert)

""" Backbone """
os.system("mkdir {output_folder}/for_rev_split_backbone".format(output_folder=output_folder))
os.chdir("{output_folder}/for_rev_split_backbone".format(output_folder=output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {cores} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_backbone/ --rf={full_reference} -g {backbone_locus} \"".format(
    timeslot=timeslot2,
    project=project,
    mem=mem,
    cores=cores,
    mail=mail,
    depend=jobid_default,
    venv=venv,
    split=split,
    output_folder=output_folder,
    full_reference=full_reference,
    backbone_locus=backbone_locus
)
jobid =subprocess.getoutput(action)
jobid_split_backbone = jobid.split()[3]
print("Submitted forward reverse split backbone\n",action,jobid_split_backbone)


""" Repeat count analysis """
""" Insert """
os.system("mkdir {output_folder}/split_insert".format(output_folder=output_folder))
os.chdir("{output_folder}/split_insert".format(output_folder=output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {cores} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_insert --rf={full_reference} -a {insert_locus}\"".format(
    timeslot=timeslot3,
    project=project,
    mem=mem,
    cores=cores,
    mail=mail,
    depend=jobid_default,
    venv=venv,
    repeat=repeat,
    output_folder=output_folder,
    full_reference=full_reference,
    insert_locus = insert_locus
)
jobid =subprocess.getoutput(action)
jobid_repeat_insert = jobid.split()[3]
print("Submitted repeat analysis insert\n",action,jobid_repeat_insert)

""" Backbone """
os.system("mkdir {output_folder}/split_backbone".format(output_folder=output_folder))
os.chdir("{output_folder}/split_backbone".format(output_folder=output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {cores} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_backbone --rf={full_reference} -a {backbone_locus}\"".format(
    timeslot=timeslot3,
    project=project,
    mem=mem,
    cores=cores,
    mail=mail,
    depend=jobid_default,
    venv=venv,
    repeat=repeat,
    output_folder=output_folder,
    full_reference=full_reference,
    backbone_locus = backbone_locus
)
jobid =subprocess.getoutput(action)
jobid_repeat_backbone = jobid.split()[3]
print("Submitted repeat analysis backbone\n",action,jobid_repeat_backbone)


"""Count alleles """
if not os.path.isdir("{output_folder}/SH/".format(output_folder=output_folder)):
    os.system("mkdir {output_folder}/SH/".format(output_folder=output_folder))

write_file=open(str(output_folder) + "/SH/Count_alleles.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/SH/Count_alleles.output\n#SBATCH -e {output_folder}/SH/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot=timeslot1,
    project=project,
    mem=mem,
    output_folder=output_folder,
    mail=mail
))
write_file.write("source {venv}".format(venv=venv))
write_file.write("cd {output_folder}".format(output_folder=output_folder))
write_file.write("{calculate}".format(calculate=calculate))
write_file.write("cd {output_folder}/for_rev_split_backbone/".format(output_folder=output_folder))
write_file.write("{calculate}".format(calculate=calculate))
write_file.write("cd {output_folder}/for_rev_split_insert/".format(output_folder=output_folder))
write_file.write("{calculate}".format(calculate=calculate))
write_file.write("cd {output_folder}/split_insert/bin_consensus/".format(output_folder=output_folder))
write_file.write("{calculate}".format(calculate=calculate))
write_file.write("rm {output_folder}/split_insert/bin_consensus_folder/".format(output_folder=output_folder))
write_file.write("cd {output_folder}/split_backbone/bin_consensus/".format(output_folder=output_folder))
write_file.write("{calculate}".format(calculate=calculate))
write_file.write("rm {output_folder}/split_backbone/bin_consensus_folder/".format(output_folder=output_folder))
write_file.close()
os.system("sbatch --depend={default},{split_insert},{split_backbone},{repeat_insert},{repeat_backbone} {output_folder}/SH/Count_alleles.sh".format(
    default=jobid_default,
    split_insert=jobid_split_insert,
    split_backbone=jobid_split_backbone,
    repeat_insert=jobid_repeat_insert,
    repeat_backbone=jobid_repeat_backbone,
    output_folder=output_folder
))
print("Submitted count alleles\n")


""" Make Structure file and plot """
write_file=open(str(output_folder) + "/SH/Make_structure.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/SH/Count_alleles.output\n#SBATCH -e {output_folder}/SH/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot = timeslot3,
    project=project,
    mem=mem,
    output_folder=output_folder,
    mail=mail,
))
write_file.write("source {venv}".format(venv=venv))
write_file.write("mkdir {output_folder}/structure".format(output_folder=output_folder))
write_file.write("cd {output_folder}/structure".format(output_folder=output_folder))
write_file.write("{structure} -i {output_folder}/bam -o {overlap} i {insert_targetinterval} > structure.txt".format(
    structure=structure,
    output_folder=output_folder,
    overlap=overlap,
    insert_targetinterval = insert_targetinterval
))
write_file.write("{rscript} {plot} structure.txt {number_reads}".format(
    rscript=rscript,
    plot=plot_dashboard,
    number_reads=number_reads
))
write_file.close()
os.system("sbatch --depend={default} {output_folder}/SH/Make_structure.sh".format(
    default=jobid_default,
    output_folder=output_folder
))
print("Submitted structure file\n")


""" Check numbers and make overview """
write_file=open(str(output_folder) + "/SH/Check.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/SH/Count_alleles.output\n#SBATCH -e {output_folder}/SH/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot = timeslot2,
    project=project,
    mem=mem,
    output_folder=output_folder,
    mail=mail,
))
write_file.write("source {venv}".format(venv=venv))
write_file.write("cd {output_folder}".format(output_folder=output_folder))
write_file.write("{check_numbers} -r {input_folder} -p {output_folder} > check_numbers.txt ".format(
    check_numbers=check_numbers,
    input_folder=input_folder,
    output_folder=output_folder,
    overlap=overlap,
    insert_targetinterval = insert_targetinterval
))
write_file.write("cd {output_folder}/bam".format(output_folder=output_folder))
write_file.write("{find_read_bam}  > bam_locations.txt ".format(
    find_read_bam=find_read_bam
))
write_file.close()
os.system("sbatch --depend={default} {output_folder}/SH/Make_structure.sh".format(
    default=jobid_default,
    output_folder=output_folder
))
print("Submitted cleanup\n")
os.chdir("{output_folder}".format(output_folder=output_folder))
