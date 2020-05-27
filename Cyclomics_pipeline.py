#! /usr/bin/env python
import sys
import os
import subprocess
import settings
import argparse

# Make into arg.parse
parser = argparse.ArgumentParser()
subparser = parser.add_subparsers()
parser_slurm = subparser.add_parser('slurm', help='Submit jobs through SLURM')
parser_slurm.add_argument('input_folder', help='Input folder')
parser_slurm.add_argument('output_folder', help='Output folder')
parser_slurm.add_argument('prefix', help='Prefix (e.g. RunID)')
parser_slurm.add_argument('insert_locus', help='Insert locus (e.g. TP53)')
parser_slurm.add_argument('backbone_locus', help='Backbone locus (e.g. BB25)')
parser_slurm.add_argument('--insert_targetinterval', default="TP53:1-27760", help='Target locus interval for structure file [default =TP53:1-27760]')
args = parser.parse_args()


"""Default consensus calling"""
action = "source {venv} && {default} -i {input_folder} -o {output_folder} -p {prefix} --rf={full_reference} --rt={target_reference}".format(
    venv=settings.venv,
    default=settings.default,
    input_folder=args.input_folder,
    output_folder=args.output_folder,
    prefix=args.prefix,
    full_reference=settings.full_ref,
    target_reference=settings.target_ref
)
os.system(action)
jobid_default = open("{output_folder}/SH/job_id_run_dagcon_consensus.sh".format(output_folder=args.output_folder), "r").readline().strip()
print("Submitted default consensus calling\n",action,jobid_default)


""" Forward and Reverse splitting """
""" Insert """
os.system("mkdir {output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
os.chdir("{output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_insert/ --rf={full_reference} -g {insert_locus} \"".format(
    timeslot=settings.SLURM_JOB_TIME_MED,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    threads=settings.THREADS,
    mail=settings.mail,
    depend=jobid_default,
    venv=settings.venv,
    split=settings.split,
    output_folder=args.output_folder,
    full_reference=settings.full_ref,
    insert_locus=args.insert_locus
)
jobid =subprocess.getoutput(action)
jobid_split_insert = jobid.split()[3]
print("Submitted forward reverse split insert\n",action,jobid_split_insert)

""" Backbone """
os.system("mkdir {output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
os.chdir("{output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_backbone/ --rf={full_reference} -g {backbone_locus} \"".format(
    timeslot=settings.SLURM_JOB_TIME_MED,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    threads=settings.THREADS,
    mail=settings.mail,
    depend=jobid_default,
    venv=settings.venv,
    split=settings.split,
    output_folder=args.output_folder,
    full_reference=settings.full_ref,
    backbone_locus=args.backbone_locus
)
jobid =subprocess.getoutput(action)
jobid_split_backbone = jobid.split()[3]
print("Submitted forward reverse split backbone\n",action,jobid_split_backbone)


""" Repeat count analysis """
""" Insert """
os.system("mkdir {output_folder}/split_insert".format(output_folder=args.output_folder))
os.chdir("{output_folder}/split_insert".format(output_folder=args.output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_insert --rf={full_reference} -a {insert_locus}\"".format(
    timeslot=settings.SLURM_JOB_TIME_HIGH,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    threads=settings.THREADS,
    mail=settings.mail,
    depend=jobid_default,
    venv=settings.venv,
    repeat=settings.repeat,
    output_folder=args.output_folder,
    full_reference=settings.full_ref,
    insert_locus = args.insert_locus
)
jobid =subprocess.getoutput(action)
jobid_repeat_insert = jobid.split()[3]
print("Submitted repeat analysis insert\n",action,jobid_repeat_insert)

""" Backbone """
os.system("mkdir {output_folder}/split_backbone".format(output_folder=args.output_folder))
os.chdir("{output_folder}/split_backbone".format(output_folder=args.output_folder))
action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --depend={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_backbone --rf={full_reference} -a {backbone_locus}\"".format(
    timeslot=settings.SLURM_JOB_TIME_HIGH,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    threads=settings.THREADS,
    mail=settings.mail,
    depend=jobid_default,
    venv=settings.venv,
    repeat=settings.repeat,
    output_folder=args.output_folder,
    full_reference=settings.full_ref,
    backbone_locus=args.backbone_locus
)
jobid =subprocess.getoutput(action)
jobid_repeat_backbone = jobid.split()[3]
print("Submitted repeat analysis backbone\n",action,jobid_repeat_backbone)


"""Count alleles """
if not os.path.isdir("{output_folder}/jobs/".format(output_folder=args.output_folder)):
    os.system("mkdir {output_folder}/jobs/".format(output_folder=args.output_folder))

write_file=open(str(args.output_folder) + "/jobs/Count_alleles.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/jobs/Count_alleles.output\n#SBATCH -e {output_folder}/jobs/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot=settings.SLURM_JOB_TIME_LOW,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    output_folder=args.output_folder,
    mail=settings.mail
))
write_file.write("source {venv}\n".format(venv=settings.venv))
write_file.write("cd {output_folder}\n".format(output_folder=args.output_folder))
write_file.write("{calculate}\n".format(calculate=settings.calculate))
write_file.write("cd {output_folder}/for_rev_split_backbone/\n".format(output_folder=args.output_folder))
write_file.write("{calculate}\n".format(calculate=settings.calculate))
write_file.write("cd {output_folder}/for_rev_split_insert/\n".format(output_folder=args.output_folder))
write_file.write("{calculate}\n".format(calculate=settings.calculate))
write_file.write("cd {output_folder}/split_insert/bin_consensus/\n".format(output_folder=args.output_folder))
write_file.write("{calculate}".format(calculate=settings.calculate))
write_file.write("rm {output_folder}/split_insert/bin_consensus_folder/\n".format(output_folder=args.output_folder))
write_file.write("cd {output_folder}/split_backbone/bin_consensus/\n".format(output_folder=args.output_folder))
write_file.write("{calculate}\n".format(calculate=settings.calculate))
write_file.write("rm {output_folder}/split_backbone/bin_consensus_folder/\n".format(output_folder=args.output_folder))
write_file.close()
os.system("sbatch --depend={default},{split_insert},{split_backbone},{repeat_insert},{repeat_backbone} {output_folder}/jobs/Count_alleles.sh".format(
    default=jobid_default,
    split_insert=jobid_split_insert,
    split_backbone=jobid_split_backbone,
    repeat_insert=jobid_repeat_insert,
    repeat_backbone=jobid_repeat_backbone,
    output_folder=args.output_folder
))
print("Submitted count alleles\n")


""" Make Structure file and plot """
write_file=open(str(args.output_folder) + "/jobs/Make_structure.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/jobs/Count_alleles.output\n#SBATCH -e {output_folder}/jobs/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot = settings.SLURM_JOB_TIME_HIGH,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    output_folder=args.output_folder,
    mail=settings.mail,
))
write_file.write("source {venv}\n".format(venv=settings.venv))
write_file.write("mkdir {output_folder}/structure\n".format(output_folder=args.output_folder))
write_file.write("cd {output_folder}/structure\n".format(output_folder=args.output_folder))
write_file.write("{structure} -i {output_folder}/bam -o {overlap} i {insert_targetinterval} > structure.txt\n".format(
    structure=settings.structure,
    output_folder=args.output_folder,
    overlap=settings.STRUCTURE_OVERLAP,
    insert_targetinterval=args.insert_targetinterval
))
write_file.write("{rscript} {plot} structure.txt {number_reads}\n".format(
    rscript=settings.rscript,
    plot=settings.plot_dashboard,
    number_reads=settings.MAX_READS_JOB_PLOT
))
write_file.close()

os.system("sbatch --depend={default} {output_folder}/jobs/Make_structure.sh".format(
    default=jobid_default,
    output_folder=args.output_folder
))
print("Submitted structure file\n")


""" Check numbers and make overview """
write_file=open(str(args.output_folder) + "/jobs/Check.sh","w")
write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/jobs/Count_alleles.output\n#SBATCH -e {output_folder}/jobs/Count_alleles.error \n#SBATCH --mail-user={mail}\n".format(
    timeslot = settings.SLURM_JOB_TIME_MED,
    project=settings.project,
    mem=settings.MAX_MEM_TARGET,
    output_folder=args.output_folder,
    mail=settings.mail,
))
write_file.write("source {venv}\n".format(venv=settings.venv))
write_file.write("cd {output_folder}\n".format(output_folder=args.output_folder))
write_file.write("{check_numbers} -r {input_folder} -p {output_folder} > check_numbers.txt\n".format(
    check_numbers=settings.check_numbers,
    input_folder=args.input_folder,
    output_folder=args.output_folder
))
write_file.write("cd {output_folder}/bam\n".format(output_folder=args.output_folder))
write_file.write("{find_read_bam}  > bam_locations.txt\n".format(
    find_read_bam=settings.find_read_bam
))
write_file.close()
os.system("sbatch --depend={default} {output_folder}/jobs/Make_structure.sh".format(
    default=jobid_default,
    output_folder=args.output_folder
))
print("Submitted cleanup\n")
os.chdir("{output_folder}".format(output_folder=args.output_folder))
