#! /usr/bin/env python
import sys
import os
import subprocess
import settings
import argparse

def slurm (args):
    """Default consensus calling"""
    action = "source {venv} && {default} -i {input_folder} -o {output_folder} -p {prefix} --rf={full_reference} --rt={target_reference}".format(
        venv=settings.venv,
        default=settings.default_slurm,
        input_folder=args.input_folder,
        output_folder=args.output_folder,
        prefix=args.prefix,
        full_reference=settings.full_ref,
        target_reference=settings.target_ref
    )
    os.system(action)
    jobid_default = open("{output_folder}/jobs/job_id_run_dagcon_consensus.sh".format(output_folder=args.output_folder), "r").readline().strip()
    print("Submitted default consensus calling\n",action,jobid_default,"\n")


    """ Forward and Reverse splitting """
    """ Insert """
    os.system("mkdir {output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
    action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --dependency={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_insert/ --rf={full_reference} -g {insert_locus} --id {sample_id}\"".format(
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
        insert_locus=args.insert_locus,
        sample_id = args.prefix
    )
    jobid = subprocess.getoutput(action)
    jobid_split_insert = jobid.split()[3]
    print("Submitted forward reverse split insert\n",action,jobid_split_insert,"\n")


    """ Backbone """
    os.system("mkdir {output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
    action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --dependency={depend} --wrap=\"source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_backbone/ --rf={full_reference} -g {backbone_locus} --id {sample_id}\"".format(
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
        backbone_locus=args.backbone_locus,
        sample_id = args.prefix
    )
    jobid = subprocess.getoutput(action)
    jobid_split_backbone = jobid.split()[3]
    print("Submitted forward reverse split backbone\n",action,jobid_split_backbone,"\n")


    """ Repeat count analysis """
    """ Insert """
    os.system("mkdir {output_folder}/split_insert".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/split_insert".format(output_folder=args.output_folder))
    action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --dependency={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_insert --rf={full_reference} -a {insert_locus} --maxfiles={maxfilecount}\"".format(
        timeslot=settings.SLURM_JOB_TIME_LONG,
        project=settings.project,
        mem=settings.MAX_MEM_TARGET,
        threads=settings.THREADS,
        mail=settings.mail,
        depend=jobid_default,
        venv=settings.venv,
        repeat=settings.repeat_slurm,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        insert_locus = args.insert_locus,
        maxfilecount=args.maxfilecount
    )
    jobid = subprocess.getoutput(action)
    jobid_repeat_insert = jobid.split()[3]
    print("Submitted repeat analysis insert\n",action,jobid_repeat_insert,"\n")


    """ Backbone """
    os.system("mkdir {output_folder}/split_backbone".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/split_backbone".format(output_folder=args.output_folder))
    action = "sbatch -t {timeslot} -A {project} --export=NONE --mem={mem}G -c {threads} --mail-type=FAIL --mail-user={mail} --dependency={depend} --wrap=\"source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_backbone --rf={full_reference} -a {backbone_locus} --maxfiles={maxfilecount}\"".format(
        timeslot=settings.SLURM_JOB_TIME_LONG,
        project=settings.project,
        mem=settings.MAX_MEM_TARGET,
        threads=settings.THREADS,
        mail=settings.mail,
        depend=jobid_default,
        venv=settings.venv,
        repeat=settings.repeat_slurm,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        backbone_locus=args.backbone_locus,
        maxfilecount=args.maxfilecount
    )
    jobid = subprocess.getoutput(action)
    jobid_repeat_backbone = jobid.split()[3]
    print("Submitted repeat analysis backbone\n",action,jobid_repeat_backbone,"\n")


    """ Make Structure file and plot """
    write_file=open(str(args.output_folder) + "/jobs/Make_structure.sh","w")
    write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/jobs/Make_structure.output\n#SBATCH -e {output_folder}/jobs/Make_structure.error \n#SBATCH --mail-user={mail}\n".format(
        timeslot = settings.SLURM_JOB_TIME_LONG,
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
        number_reads=args.structure_plot_max
    ))
    write_file.close()

    os.system("sbatch --dependency={default} {output_folder}/jobs/Make_structure.sh".format(
        default=jobid_default,
        output_folder=args.output_folder
    ))
    print("Submitted structure file\n")


    """ Check numbers and make overview of reads in targeted BAMs tars """
    write_file=open(str(args.output_folder) + "/jobs/Check.sh","w")
    write_file.write("#!/bin/bash\n#SBATCH -t {timeslot} \n#SBATCH --account={project} \n#SBATCH --mem={mem}G \n#SBATCH --export=NONE\n#SBATCH -o {output_folder}/jobs/Check.output\n#SBATCH -e {output_folder}/jobs/Check.error \n#SBATCH --mail-user={mail}\n".format(
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
    os.system("sbatch --dependency={default} {output_folder}/jobs/Check.sh".format(
        default=jobid_default,
        output_folder=args.output_folder
    ))
    print("Submitted cleanup\n")
    os.chdir("{output_folder}".format(output_folder=args.output_folder))




def nocluster (args):

    """Default consensus calling"""
    action = "source {venv} && {default} -i {input_folder} -o {output_folder} -p {prefix} --rf={full_reference} --rt={target_reference}".format(
        venv=settings.venv,
        default=settings.default_nocluster,
        input_folder=args.input_folder,
        output_folder=args.output_folder,
        prefix=args.prefix,
        full_reference=settings.full_ref,
        target_reference=settings.target_ref
    )
    print("Default consensus calling:",action)
    os.system(action)
    print("Finished default consensus calling\n")

    """ Forward and Reverse splitting """
    """ Insert """
    os.system("mkdir -p {output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/for_rev_split_insert".format(output_folder=args.output_folder))
    action = "source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_insert/ --rf={full_reference} -g {insert_locus} --id {sample_id}".format(
        venv=settings.venv,
        split=settings.split,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        insert_locus=args.insert_locus,
        sample_id = args.prefix
    )
    print("Forward and Reverse splitting, Insert:",action)
    os.system(action)
    print("Finshed forward reverse split insert\n")


    """ Backbone """
    os.system("mkdir {output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/for_rev_split_backbone".format(output_folder=args.output_folder))
    action = "source {venv} && {split} -i {output_folder}/bam -o {output_folder}/for_rev_split_backbone/ --rf={full_reference} -g {backbone_locus} --id {sample_id}".format(
        venv=settings.venv,
        split=settings.split,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        backbone_locus=args.backbone_locus,
        sample_id = args.prefix
    )
    print("Forward and Reverse splitting, Backbone:", action)
    os.system(action)
    print("Finished forward reverse split backbone\n")


    """ Repeat count analysis """
    """ Insert """
    os.system("mkdir {output_folder}/split_insert".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/split_insert".format(output_folder=args.output_folder))
    action = "source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_insert --rf={full_reference} -a {insert_locus}".format(
        venv=settings.venv,
        repeat=settings.repeat_nocluster,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        insert_locus = args.insert_locus
    )
    print("Repeat count analysis, Insert:", action)
    os.system(action)
    print("Finshed repeat analysis insert\n")


    """ Backbone """
    os.system("mkdir {output_folder}/split_backbone".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/split_backbone".format(output_folder=args.output_folder))
    action = "source {venv} && {repeat} --ib {output_folder}/bam/ --im {output_folder}/m5/ -o {output_folder}/split_backbone --rf={full_reference} -a {backbone_locus}".format(
        venv=settings.venv,
        repeat=settings.repeat_nocluster,
        output_folder=args.output_folder,
        full_reference=settings.full_ref,
        backbone_locus=args.backbone_locus
    )
    print("Repeat count analysis, Backbone:", action)
    os.system(action)
    print("Finshed repeat analysis backbone\n")


    """ Make Structure file and plot """
    os.system("mkdir {output_folder}/structure\n".format(output_folder=args.output_folder))
    os.chdir("{output_folder}/structure".format(output_folder=args.output_folder))
    action = "source {venv} && {structure} -i {output_folder}/bam -o {overlap} i {insert_targetinterval} > structure.txt".format(
        venv=settings.venv,
        structure=settings.structure,
        output_folder=args.output_folder,
        overlap=settings.STRUCTURE_OVERLAP,
        insert_targetinterval=args.insert_targetinterval
    )
    print("Make Structure file and plot:", action)
    os.system(action)

    action = ("{rscript} {plot} structure.txt {number_reads}\n".format(
        rscript=settings.rscript,
        plot=settings.plot_dashboard,
        number_reads=args.structure_plot_max
    ))
    print(action)
    os.system(action)
    print("Finished structure file\n")


    """ Check numbers and make overview of reads in targeted BAMs tars """
    os.chdir("{output_folder}".format(output_folder=args.output_folder))
    action = "source {venv} && {check_numbers} -r {input_folder} -p {output_folder} > check_numbers.txt\n".format(
        venv=settings.venv,
        check_numbers=settings.check_numbers,
        input_folder=args.input_folder,
        output_folder=args.output_folder
    )
    print(action)
    os.system(action)

    os.chdir("{output_folder}/bam".format(output_folder=args.output_folder))
    action = "source {venv} && {find_read_bam}  > bam_locations.txt\n".format(
        venv=settings.venv,
        find_read_bam=settings.find_read_bam
    )
    print(action)
    os.system(action)

    print("Finished cleanup\n")
    os.chdir("{output_folder}".format(output_folder=args.output_folder))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    parser_slurm = subparser.add_parser('slurm', help='submit parallel jobs with SLURM')
    parser_slurm.add_argument('input_folder', help='Input folder')
    parser_slurm.add_argument('output_folder', help='Output folder')
    parser_slurm.add_argument('prefix', help='Prefix (e.g. RunID/SampleID)')
    parser_slurm.add_argument('insert_locus', help='Insert locus (e.g. TP53)')
    parser_slurm.add_argument('backbone_locus', help='Backbone locus (e.g. BB25)')
    parser_slurm.add_argument('--insert_targetinterval', default="TP53:1-27760", help='Target locus interval for structure file [default =TP53:1-27760]')
    parser_slurm.add_argument('--structure_plot_max', default=settings.MAX_READS_JOB_PLOT, help='Maximum reads in structure file plotted [default MAX_READS_JOB_PLOT in settings.py ]')
    parser_slurm.add_argument('--maxfilecount', default=settings.MAX_FILE_COUNT, help='Maximum files used in repeat counting [default MAX_FILE_COUNT in settings.py]')
    parser_slurm.set_defaults(func = slurm)

    parser_nocluster = subparser.add_parser('nocluster', help='do not use parallel jobs (commandline only)')
    parser_nocluster.add_argument('input_folder', help='Input folder')
    parser_nocluster.add_argument('output_folder', help='Output folder')
    parser_nocluster.add_argument('prefix', help='Prefix (e.g. RunID/SampleID)')
    parser_nocluster.add_argument('insert_locus', help='Insert locus (e.g. TP53)')
    parser_nocluster.add_argument('backbone_locus', help='Backbone locus (e.g. BB25)')
    parser_nocluster.add_argument('--insert_targetinterval', default="TP53:1-27760", help='Target locus interval for structure file [default =TP53:1-27760]')
    parser_nocluster.add_argument('--structure_plot_max', default=settings.MAX_READS_JOB_PLOT, help='Maximum reads in structure file plotted [default MAX_READS_JOB_PLOT in settings.py ]')
    parser_nocluster.set_defaults(func = nocluster)

    args = parser.parse_args()
    # Temp fix for relative paths but probably breaks full paths
    # Python lacks a proper realpath, MacOS completely lacks a realpath.
    # "Think different", yo.
    args.input_folder  = os.path.abspath(args.input_folder)
    args.output_folder = os.path.abspath(args.output_folder)
    args.func(args)
