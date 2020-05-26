
mail = 'm.elferink@umcutrecht.nl'
project = 'compgen'


venv = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/venv/bin/activate"
bwa = "/hpc/local/CentOS7/cog_bioinf/bwa-0.7.17/bwa"
sambamba = "/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
bam2m5 = "/hpc/compgen/tools/bam2m5/bam2m5.py"
pbdagcon = "/hpc/compgen/tools/pbdagcon/src/cpp/pbdagcon"
last = "/hpc/compgen/tools/last-921/"  #Note that this is the folder. Folder should contain ./src/lastal, ./src/last-split, ./scripts/maf-convert, and ./last_params

last_param = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/last_params"
#last_param = "/hpc/compgen/tools/development/DEV_Cyclomics_consensus_pipeline/Test_Last/flongle_DER4803.param"

full_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta"
target_ref = "/hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/BRAF_TP53_EGFR_BB_pjet.fasta"

DAGCON_MIN_COV = 10
SLURM_JOB_TIME = "4:00:00"
SLURM_PARALLEL_JOBS = 500 
MAX_READS_JOB = 1000 
MAX_MEM_TARGET = 10
MAX_MEM_FULL  = 32 
THREADS = 2 
MIN_CONS_LEN = 35
TRIM = 0
 
#split_forward_reverse_reads.py
SPLIT_MIN_PERC = 50 

#bin_on_repeat_count.py 
INSERT = "TP53"
PBDAGCON_PARAM = " -c 1 -t 0 -j 2 "

#make_structure.py 
STRUCTURE_OVERLAP = 100
STRUCTURE_FLANK = 30
STRUCTURE_MAD = 100
STRUCTURE_INTARGET = "17:7565097-7590856"

#filter_sam.py
samtools = "/hpc/local/CentOS7/cog/software/samtools-1.2/samtools"
