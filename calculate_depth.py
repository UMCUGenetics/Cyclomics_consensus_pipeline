#! /usr/bin/env python
import sys, os,re, commands
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-k", dest="target_key", metavar="[PATH]", help="target_key [e.g. TP53]")
        group.add_option("-v", dest="target_value", metavar="[PATH]", help="value_key [e.g. 17:7565097-7590856]")
        group.add_option("-r", default="yes", dest="cosmic", help="calculate COSMIC position TP53 yes/no [default= yes]")
        group.add_option("-b", dest="blacklist", help="blacklist of reads that should be excluded for basecalling")
        group.add_option("-o", dest="old", help="calculates stats old reference contigs [default= off]")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

dic={}
if opt.target_key and opt.target_value:
    dic[opt.target_key]=opt.target_value
else:
    if opt.old:
        dic ={"exon12":"17:7577000-7577180","TP53":"17:7565097-7590856","BB200_1":"BB200_1:1-243","BB200_2":"BB200_2:1-244","BB200_3":"BB200_3:1-244","BB200_4n":"BB200_4n:1-247","BB200_5":"BB200_5:1-243","pJet":"pJet:1-2974"}
    else:
        dic={"exon12":"17:7577000-7577180","BB22":"BB22:1-248","BB24":"BB24:1-248","BB25":"BB25:1-248","BBCR":"BB22:1-248","TP53":"17:7565097-7590856","PJET":"PJET:1-2974","EGFR":"7:55081714-55329313"}


sambamba="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
cwd = os.getcwd()
cosmic= "/hpc/compgen/tools/Cyclomics_consensus_pipeline/data_files/COSMIC_mutations.bed"

for f in os.listdir(cwd):
    if f.endswith((".sorted.bam")):
        if opt.blacklist:
            f="{0}_bl.bam".format(str(f)[0:-4])
 
        for item in dic:
            if opt.blacklist:
                new_file=str(f).replace('.sorted_bl.bam','_sambamba_output_bl_'+str(item)+'.txt')
            else:
                new_file=str(f).replace('.sorted.bam','_sambamba_output_'+str(item)+'.txt')

            action = "{0} depth base -L {1} --min-coverage=0 {2} > {3} ".format(
                sambamba,
                str(dic[item]),
                str(f),
                new_file
            )
            os.system(action)
        if opt.cosmic=="yes":
            if opt.blacklist:
                new_file=str(f).replace('.sorted_bl.bam','_sambamba_output_bl_cosmic.txt')
            else:
                new_file=str(f).replace('.sorted.bam','_sambamba_output_cosmic.txt')

            action = "{0} depth base -L {1} --min-coverage=0 {2} > {3} ".format(
                sambamba,
                str(cosmic),
                str(f),
                new_file
            )
            os.system(action)
