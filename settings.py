# User specific details
mail = 'm.elferink@umcutrecht.nl'
project = 'compgen'


# Paths to scripts
repo_path = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/"
venv = str(repo_path)+ "venv/bin/activate"
default = str(repo_path) + "/run_dagcon_consensus.py"
split = str(repo_path) + "/split_forward_reverse_reads.py"
repeat = str(repo_path) + "/bin_on_repeat_count.py"
calculate = str(repo_path) + "/calculate_depth.py"
plot_dashboard = str(repo_path) + "/plot_Dashboard.R"
check_numbers = str(repo_path) + "/check_numbers.py"
find_read_bam = str(repo_path) + "/find_read_bam.py"
structure = str(repo_path) + "/make_structure.py"

# Paths to third party tools
bwa = "/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa"
sambamba = "/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
bam2m5 = "/hpc/compgen/tools/bam2m5/bam2m5.py"
pbdagcon = "/hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon"
lastal = "/hpc/compgen/tools/last-921/src/lastal"
lastsplit = "/hpc/compgen/tools/last-921/src/last-split"
lastparam = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/data_files/last_params"
mafconvert = "/hpc/compgen/tools/last-921/scripts/maf-convert"
rscript= "/hpc/local/CentOS7/common/lang/R/3.2.2/bin/Rscript"


# Reference Genome files 
full_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta"
target_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta"


#  General settings


SLURM_JOB_TIME_LOW = "2:00:00"
SLURM_JOB_TIME_MED = "4:00:00"
SLURM_JOB_TIME_HIGH = "196:00:00"


MAX_MEM_TARGET = 10
MAX_MEM_FULL  = 32 
THREADS = 2 
MIN_CONS_LEN = 35
TRIM = 0
BWA_MEM = "{bwa} mem -t 2 -c 100 -M -R".format(bwa=bwa)
MAX_READS_JOB_PLOT = 100000




#run_dagcon_consensus.py specific settings
DAGCON_MIN_COV = 10
SLURM_PARALLEL_JOBS = 500
MAX_READS_JOB = 1000

#split_forward_reverse_reads.py specific settings
SPLIT_MIN_PERC = 50 
STATES = ["forward","reverse"]

#bin_on_repeat_count.py specific settings
INSERT = "TP53"
PBDAGCON_PARAM = " -c 1 -t 0 -j 2 "

#make_structure.py specific settings
STRUCTURE_OVERLAP = 100
STRUCTURE_FLANK = 30
STRUCTURE_MAD = 100
STRUCTURE_INTARGET = "17:7565097-7590856"

#filter_sam.py specific settings
samtools = "/hpc/local/CentOS7/cog/software/samtools-1.2/samtools"
