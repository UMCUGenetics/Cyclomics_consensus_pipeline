# Cyclomics_consensus_pipeline
Collection of scripts to process CyclomicsSeq (nanopore) data.

## Locally install these tools (version are tested versions
bwa	v0.7.17	https://github.com/lh3/bwa \
sambamba	v0.6.5	https://github.com/biod/sambamba \
bam2m5	commit=0ef1a930b6a0426c55e8de950bf1ac22eef61bdf	https://github.com/sein-tao/bam2m5 \
pbdagcon	tag=p4-mainline-127-g3c382f2	commit=3c382f2673fbf3c5305f5323188e790dc396ac9d	https://github.com/PacificBiosciences/pbdagcon \
last-921	v921	http://last.cbrc.jp/ \
R/Rscript	v3.2.2	https://www.r-project.org/ 

## Make virtual python environment
_ _(Tested with python v3.6.1)_ _ \
virtualenv -p python3 venv \
source venv/bin/activate \
easy_install pip \
pip install -r requirements.txt 

## Index full and targeted reference genomes
Picard	(tested=v1.141) \
lastdb	(v921) \
bwa	(v0.7.17)
 
## Run Scripts
__Always load virtualenv before running scripts__ 
```bash
source venv/bin/activate
```

1) Cyclomics_pipeline.py
Wrapper script that includes scripts 2-9 as described below. \
    slurm            submit parallel jobs with SLURM. \
    nocluster        do not use parallel jobs (commandline only) \

Usage:
```bash
python Cyclomics_pipeline.py {slurm/nocluster} {raw_data folder with fastq files} {output folder} {prefix (eg run or sampleID)} {insert locus (e.g. TP53)} {backbone locus (e.g. BB25)}
```
_ _optional:_ _  \
    for either slurm or nocluster: \
        --insert_targetinterval   	structure file: define what is considered in-target \
        --structure_plot_max 		structure plot: maximum number of reads included in the plot. Note that this should be less than the number of reads in the run. \
    for slurm only: \
        --maxfilecount 			bin_on_repeat_count.py: maximum number of files used in bin repeat. This might be helpful with large runs that will take a very long time if all data is used. Number of file * number MAX_READS_JOB = max number of reads.

	
# Script that are used in Cyclomics_pipeline.py, but can also be manually runned
2) run_dagcon_consensus.py / run_dagcon_consensus_nocluster.py \
This is the main script to process the Cyclomics nanopore data.  
Nanopore reads will be LAST-split to a targeted reference genome including the backbone and insert(gene) sequences.  
PBDAGCON is used to create consensus sequences for backbone and insert repeats (with a minimum threshold op repeats needed). 
These consensus reads will be mapped to the full reference genome using BWA.

3) split_forward_reverse_reads.py \
Script that will divide reads into forward or reverse sequenced reads (for backbone or insert).
This analysis can be performed to determine if one strand performs better than the other for specific basepositions.

4) bin_on_repeat_count.py / bin_on_repeat_count_nocluster.py \
Script that will divide reads into specific repeat bins (for backbone or insert).
This analysis can be performed to see the influence of increased number of repeats.

5) calculate_depth.py \
Script that will count the alleles for a specific target based on BAM files

6) make_structure.py \
Script that will use the targeted BAM to determine the original structure of the full read with regards to insert and backbone.

7) check_numbers.py \
Script to compare reads numbers in raw data, bam, m5, consensus folders. To make sure jobs were killed or not finished.

8) find_read_bam.py \
Script to make a single file in which readID are linked to the tar-ball. Useful to quickly extract specific reads if needed.

9) plot_Dashboard.R \
Rscript to plot statistics from the make_structure.py output file.

#Additional scripts (not used in Cyclomics_pipeline.py) 

10) run_dagcon_consensus_nocluster.py \
Like run_dagcon_consensus.py, but without the use of a scheduler. Note that this might take a long time to process the data.

11) filter_sam.py \
Script to exclude specific reads in a BAM file.
