# Cyclomics_consensus_pipeline
Collection of scripts to process Cyclomics nanopore data.

## Locally install these tools:
bwa	v0.7.17	https://github.com/lh3/bwa
sambamba	v0.6.5 https://github.com/biod/sambamba
bam2m5	commit=0ef1a930b6a0426c55e8de950bf1ac22eef61bdf		https://github.com/sein-tao/bam2m5
pbdagcon	tag=p4-mainline-127-g3c382f2, commit=3c382f2673fbf3c5305f5323188e790dc396ac9d	https://github.com/PacificBiosciences/pbdagcon
last-921	v921	http://last.cbrc.jp/

## Make virtual python enviroment
(Tested with python v3.6.1)
virtualenv -p python3 venv
source venv/bin/activate
easy_install pip
pip install -r requirements.txt


## Index full and targeted reference genomes
Picard (tested=v1.141)
lastdb (v921)
bwa (v0.7.17)
 
## Run Scripts
# Always load virtualenv before running scripts
```bash
source venv/bin/activate
```

0) Consensus_calling_wrapper.py
Wrapper script to start Scripts 1-8 as decribed below. 
This script uses Slurm as scheduler

1) run_dagcon_consensus.py
This is the main script to process the Cyclomics nanopore data. 
Nanopore reads will be LAST-split to a targeted reference genome including the backbone and insert(gene) sequences.
PBDAGCON is used to create consensus sequences for backbone and insert repeats (with a minimum threshold op repeats needed). 
These consensus reads will be mapped to the full reference genome using BWA.

``` bash
python run_dagcon_consensus.py -i {input folder with nanopore fastq} -o {output folder} -p {SampleID included in BAM file naming} --rf={Full reference genome (fasta)} --rt={Target reference genome (fasta)} 
```

2) calculate_depth.py

3) make_structure.py

4) split_forward_reverse_reads.py

5) check_numbers.py

6) bin_on_repeat_count.py

7) find_read_bam.py

8) plot_Dashboard.R
Rscript to plot statistics from the make_structure.py output file.





Additional scripts
8) run_dagcon_consensus_nocluster.py
Like run_dagcon_consensus.py, but without the use of a scheduler. Note that this might take a long time to process the data.

9) filter_sam.py
Script to exclude specific reads in a BAM file.
