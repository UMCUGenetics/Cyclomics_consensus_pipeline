# Cyclomics_consensus_pipeline
Collection of scripts to process CyclomicsSeq (nanopore) data.

## Locally install these tools (version is tested version)
git - Also use `git clone` to get this repo instead of a zip \
bwa	v0.7.17	https://github.com/lh3/bwa \
sambamba	v0.6.5	https://github.com/biod/sambamba \
samtools 1.2/ \
bam2m5	commit=0ef1a930b6a0426c55e8de950bf1ac22eef61bdf	https://github.com/sein-tao/bam2m5 \
pbdagcon	tag=p4-mainline-127-g3c382f2	commit=3c382f2673fbf3c5305f5323188e790dc396ac9d	https://github.com/PacificBiosciences/pbdagcon - See below for a guide to compile this\
last-921	v921	http://last.cbrc.jp/ \
R/Rscript	v3.2.2	https://www.r-project.org/ \
Picard	(tested=v1.141) \
Java (version=1.8.0_121) \
Python v2.x, x>=6, tested v2.7 - for last \
Python v3.6 - See next paragraph

## Required R libraries (version is tested version)
ggplot2 (3.2.1)
gridExtra (2.3)
grid (3.2.2)
stringr (1.4.0)


## Make virtual python environment
_(Tested with python v3.6.1)_
```bash
virtualenv -p python3 venv
source venv/bin/activate
easy_install pip
pip install -r requirements.txt
```

Conda version:
```bash
./conda create -n cyclocons python=3.6
source ./activate cyclocons
pip install -r requirements.txt
```

## Prepare reference genome fasta

### Targeted Genome
The targeted genome comprises the sequence(s) of the gene(s) or locus of interest (e.g. TP53) and the sequence(s) of the relevant backbone(s).
This reference genome needs to be in the FASTA format, each gene or backbone added as a seperate contig.

### Full Genome
The full reference genome comprises the full target genome (e.g. hg19 or hg38) and the sequence(s) of the relevant backbone(s).
Do not add the the gene(s) of interest if the sequence is already present in the full reference genome.
This reference genome needs to be in the FASTA format, each backbone sequence added as a seperate contig.

## Index full and targeted reference genomes
Assuming your ref genome path and filename without extension is `${ref_genome}` (example: `${ref_genome}.fasta`), these steps should setup the indexes and such for tools used:
### Picard
```bash
java -jar ${path_to_picard}/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=${ref_genome}.fasta \
  OUTPUT=${ref_genome}.dict
```
### lastdb
The `lastdb` script is part part of `last` and can thus be found in the same folder.
```bash
${path_to_last}/lastdb ${ref_genome} ${ref_genome}.fasta
```

Be sure the reference genome basename and path is correct, a mistake in the naming here will give this error later on:
`lastal: can't open file: [...].prj`

### bwa
```bash
bwa index ${ref_genome}.fasta
```

## Setting.py contains paths and parameters used in the scripts.
Change settings.py according to your own system/setup. See settings.py.
Make sure the `venv` variable is an executable command that starts the python environment prepared earlier.
The `mafconvert` variable needs to combine a `python2` call and the path to the script, be sure there's a space at the end of `python2 `.


## Run Scripts
__Always load virtualenv before running scripts__
```bash
source venv/bin/activate
```

1) cyclomics_pipeline.py
Wrapper script that includes scripts 2-9 as described below. \
    slurm            submit parallel jobs with SLURM. \
    nocluster        do not use parallel jobs (commandline only)

Usage:
```bash
python cyclomics_pipeline.py {slurm/nocluster} {raw_data folder with fastq files} {output folder} {prefix (eg run or sampleID)} {insert locus (e.g. TP53)} {backbone locus (e.g. BB25)}
```
_optional:_  \
    for either slurm or nocluster: \
        --insert_targetinterval   	structure file: define what is considered in-target \
        --structure_plot_max 		structure plot: maximum number of reads included in the plot. Note that this should be less than the number of reads in the run. \
    for slurm only: \
        --maxfilecount 			bin_on_repeat_count.py: maximum number of files used in bin repeat. This might be helpful with large runs that will take a very long time if all data is used. Number of file * number MAX_READS_JOB = max number of reads.

Note: use full paths for the input/raw_data and output folder.\
Furthermore, fastq ór fastq.gz searching will be recursively in the input/raw_data folder.





__Script that are used in cyclomics_pipeline.py, but can also be manually runned__

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

#Additional scripts (not used in cyclomics_pipeline.py)

10) run_dagcon_consensus_nocluster.py \
Like run_dagcon_consensus.py, but without the use of a scheduler. Note that this might take a long time to process the data.

11) filter_sam.py \
Script to exclude specific reads in a BAM file.


# Fixes for external tools
## pbdagcon
To make pbdagcon run these days (on a Mac), there are several issues that need to be resolved. Most of these steps are simply because the code wasn't updated for a long time, and things have changed since. These steps provide solutions for issues we ran into, several may be applicable to installing the pipeline on another Unix based system.

First of all, be sure to use `git clone` to obtain pbdagcon rather than downloading a zipped version, so it can pick up the submodules it links to.
```bash
git clone https://github.com/PacificBiosciences/pbdagcon.git
```

### blasr
On MacOS 10.12 and later these errors are likely to pop up when building pbdagcon:
```
In file included from MappingMetrics.cpp:3:
./MappingMetrics.hpp:16:5: error: redefinition of enumerator '_CLOCK_REALTIME'
    CLOCK_REALTIME,
    ^
/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/time.h:154:24: note: expanded from macro 'CLOCK_REALTIME'
#define CLOCK_REALTIME _CLOCK_REALTIME
                       ^
/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/time.h:153:1: note: previous definition is here
_CLOCK_REALTIME __CLOCK_AVAILABILITY = 0,
^
```
Basically a set of definitions that existed in Linux didn't exist in MacOS until this update. After that update, the custom definitions in this package are picking a fight with the newly provided ones.

As a fix, find and open this file:
```
./blasr_libcpp/alignment/MappingMetrics.cpp
```

And around line 15, comment out the typedef enum { ... } part, making it look like this:
```cpp
/**
typedef enum {
    CLOCK_REALTIME,
    CLOCK_MONOTONIC,
    CLOCK_PROCESS_CPUTIME_ID,
    CLOCK_THREAD_CPUTIME_ID
} clockid_t;**/
```

Solution based on this thread:
https://github.com/tpaviot/oce/issues/643


### Boost
In pbdagcon, configure.py is hardcoded to download a boost headers package from a dead dropbox link.
There's information on how to get the necessary files here:
https://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html

Simplifying the process:
1. Download the whole boost_1_58_0.tar.bz2 file linked in `1 Get Boost`.
2. Unpack and put the boost subfolder (the one that directly contains files like `align.hpp`) and put it in `/pbdagcon/src/cpp/third-party/`, making the structure like so: `/pbdagcon/src/cpp/third-party/boost/align.hpp`.
3. Update configure.py to refer to this folder instead:
```python
     uri = 'https://www.dropbox.com/s/g22iayi83p5gbbq/boost_1_58_0-headersonly.tbz2?dl=0'
-    hdir = os.path.join(build_dir, 'src', 'cpp', 'third-party', 'boost_1_58_0-headersonly')
+    hdir = os.path.join(build_dir, 'src', 'cpp', 'third-party', 'boost')#_1_58_0-headersonly')
     if not os.path.isdir(hdir):
```


### DAZZLER
Hoping that a latter version of pbdagcon would fix some things we tested (only) the latest commit (`c14c422e609a914f0139f7222202ac1bce7e3ef1`) on MacOS. This version appears to assume the DAZZLER commit it downloads is a more recent version than it actually picks up and then refers to the wrong struct definition.
This gives the following error upon compiling:
```
In file included from DazAlnProvider.cpp:10:
./DazAlnProvider.hpp:97:12: error: unknown type name 'DAZZ_DB'
    Target(DAZZ_DB& db, int tspace, int small);
           ^
```

Rather than updating the linked repo, we can just make new name refer to the old one. Looks like pbdagcon doesn't expect any altered behaviour from DAZZLER (based on commits around that time) so this is probably the safe option. Open `DazAlnProvider.hpp`, and right after this block:
```c
// Dazzler headers
extern "C" {
#include "DB.h"
#include "align.h"
}
```
Insert this line, it should become about line 18:
```c
#define DAZZ_DB HITS_DB
```

After this you should be able to follow the Build/Compile steps described in the pbdagcon readme.

## last-921
Finding this particular version requires some searching on their website, this link takes you to the right commit straight away:
http://last.cbrc.jp/last/index.cgi/rev/8430dbf5abe7
Download links for zip formats can be found near the top of the page.
