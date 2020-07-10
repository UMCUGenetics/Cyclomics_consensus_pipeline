# User specific details
mail = 'm.elferink@umcutrecht.nl'	#Email adress
project = 'compgen'			#Project name for Slurm cluster

# Paths to scripts
repo_path = "/hpc/compgen/tools/Cyclomics_consensus_pipeline/"
venv = str(repo_path)+ "venv/bin/activate"
# Don't touch, should always be the same:
default_slurm = repo_path + "run_dagcon_consensus.py"
repeat_slurm = repo_path + "bin_on_repeat_count.py"
default_nocluster = repo_path + "run_dagcon_consensus_nocluster.py"
repeat_nocluster = repo_path + "bin_on_repeat_count_nocluster.py"
split = repo_path + "split_forward_reverse_reads.py"
calculate = repo_path + "calculate_depth.py"
plot_dashboard = repo_path + "plot_Dashboard.R"
check_numbers = repo_path + "check_numbers.py"
find_read_bam = repo_path + "find_read_bam.py"
structure = repo_path + "make_structure.py"
cosmic = repo_path + "data_files/COSMIC_mutations.bed"
lastparam = repo_path + "data_files/last_params"


# Paths to tools
bwa = "/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa"
sambamba = "/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
bam2m5 = "/hpc/compgen/tools/bam2m5/bam2m5.py"
pbdagcon = "/hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon"
lastal = "/hpc/compgen/tools/last-921/src/lastal"
lastsplit = "/hpc/compgen/tools/last-921/src/last-split"
mafconvert = "/hpc/compgen/tools/last-921/scripts/maf-convert"
rscript= "/hpc/local/CentOS7/common/lang/R/3.2.2/bin/Rscript"
samtools = "/hpc/local/CentOS7/cog/software/samtools-1.2/samtools"

# Reference Genome files 
full_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta"
target_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta"


#  General settings
SLURM_JOB_TIME_SHORT = "4:00:00"        # Job time SLURM
SLURM_JOB_TIME_MED = "144:00:00"        # Job time SLURM
SLURM_JOB_TIME_LONG = "196:00:00"       # Job time SLURM

MAX_MEM_TARGET = 12			# Memory SLURM (low for targeted mapping) 
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

