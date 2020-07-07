import os

# User specific details
mail = 'm.elferink@umcutrecht.nl'	#Email adress
project = 'compgen'			#Project name for Slurm cluster

# Paths to scripts
repo_path = "/Users/rstraver/Workspace/Cyclomics_consensus_pipeline/"
venv = "/Users/rstraver/Workspace/conda_tombo/bin/activate cyclocons"
# Don't touch, should always be the same:
default_slurm = str(repo_path) + "run_dagcon_consensus.py"
repeat_slurm = str(repo_path) + "bin_on_repeat_count.py"
default_nocluster = str(repo_path) + "run_dagcon_consensus_nocluster.py"
repeat_nocluster = str(repo_path) + "bin_on_repeat_count_nocluster.py"
split = str(repo_path) + "split_forward_reverse_reads.py"
calculate = str(repo_path) + "calculate_depth.py"
plot_dashboard = str(repo_path) + "plot_Dashboard.R"
check_numbers = str(repo_path) + "check_numbers.py"
find_read_bam = str(repo_path) + "find_read_bam.py"
structure = str(repo_path) + "make_structure.py"

# Paths to tools
tool_path = repo_path + "tools/"
bwa = tool_path + "bwa"
sambamba = tool_path + "sambamba_v0.6.5"
bam2m5 = tool_path + "bam2m5-master/bam2m5.py"
pbdagcon = tool_path + "pbdagcon"
lastal = tool_path + "last-921/bin/lastal"
lastsplit = tool_path + "last-921/bin/last-split"
mafconvert = 'python2 ' + tool_path + "last-921/bin/maf-convert"
rscript= "Rscript"
samtools = tool_path + "samtools"

lastparam = repo_path + "data_files/last_params" #whut?
cosmic = repo_path + "data_files/COSMIC_mutations.bed"

# Reference Genome files
full_ref = repo_path + "/data_files/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta"
target_ref = repo_path + "/data_files/testdata/BRAF_TP53_EGFR_BB_pjet.fasta"

#  General settings
SLURM_JOB_TIME_SHORT = "2:00:00"	# Job time SLURM
SLURM_JOB_TIME_MED = "8:00:00"		# Job time SLURM
SLURM_JOB_TIME_LONG = "196:00:00"	# Job time SLURM
MAX_MEM_TARGET = 10			# Memory SLURM (low for targeted mapping)
MAX_MEM_FULL  = 32 			# Memory SLURM (higher for full mapping)
THREADS = 2 				# Threads to be used in several jobs
MIN_CONS_LEN = 35			# Minimum consensus lenght for pbdagcon
TRIM = 0				# Trim x basepairs from consensus pbdagcon
BWA_MEM = "{bwa} mem -t 2 -c 100 -M -R".format(bwa=bwa)	# settings for BWA-mem
MAX_READS_JOB_PLOT = 50000  		# Maximum reads in strucure plot. Note that this number must be lower tham the total reads in the analyses.

#run_dagcon_consensus.py specific settings
DAGCON_MIN_COV = 10			# Minimum coverage needed for consensus calling with pbdagcon
SLURM_PARALLEL_JOBS = 200		# Maximum number of parallel jobs
MAX_READS_JOB = 1000			# Number of reads per jobs

#split_forward_reverse_reads.py specific settings
SPLIT_MIN_PERC = 50 			# Minimal percentage to determine either forward or reverse mapped insert/backbone (50% = majority vote)
STATES = ["forward","reverse"]		# States of mapping.

#bin_on_repeat_count.py specific settings
INSERT = "TP53"				# Dedault locus to be considered in-target
PBDAGCON_PARAM = " -c 1 -t 0 -j 2 "	# pbdagcon settings
MAX_FILE_COUNT = 0			# Maximum number of files used in bin repeat. This might be helpful with large runs that will take a very long time if all data is used. default 0 = off

#make_structure.py specific settings
STRUCTURE_OVERLAP = 100			# Overlap in bp to determine if insert fragements are from the same molecule/amplicon
STRUCTURE_FLANK = 30			# Flank to determine unmapped regions to either BB or I
STRUCTURE_MAD = 100			# If MAD score for insert startsite is more than threshold, report threshold
STRUCTURE_INTARGET = "17:7565097-7590856"	# default repeat loci considered insert target
