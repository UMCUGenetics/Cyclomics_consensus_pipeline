#! /usr/bin/env python
import sys, os,re, commands

sambamba="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
cwd = os.getcwd()
cosmic= "/hpc/cog_bioinf/ridder/tools/Cyclomics_consensus_pipeline/data_files/COSMIC_mutations.bed"
exon="17:7577000-7577180"
bb="BB200_4n:1-248"
pjet="pJet:1-2974"

for f in os.listdir(cwd):
    if f.endswith((".sorted.bam")):
        action1= str(sambamba)+" depth base -L "+str(cosmic)+ " --min-coverage=0 " + str(f) + " > "+ str(f).replace('.sorted.bam','_sambamba_output_cosmic.txt') 
        action2= str(sambamba)+" depth base -L "+str(exon)+ " --min-coverage=0 " + str(f) + " > "+ str(f).replace('.sorted.bam','_sambamba_output_Exon12.txt')
        action3= str(sambamba)+" depth base -L "+str(bb)+ " --min-coverage=0 " + str(f) + " > "+ str(f).replace('.sorted.bam','_sambamba_output_BB200_4n.txt')
        action4= str(sambamba)+" depth base -L "+str(pjet)+ " --min-coverage=0 " + str(f) + " > "+ str(f).replace('.sorted.bam','_sambamba_output_pJet.txt')

        os.system(action1)
        os.system(action2)
        os.system(action3)
        os.system(action4)




