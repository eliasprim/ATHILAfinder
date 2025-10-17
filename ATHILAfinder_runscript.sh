#!/bin/bash

# Unzip the zipped assembly file

# gzip -d Col-CEN_v1.2.fasta.gz

# Pre-process assembly file

python scripts/treat_fasta.py $1 $2




# "Assembly: $1";

# "PBS signature: $2";

# "PPT signature: $3";

# "Searching Regions (bedtools window): $4";

# "Minimum Internal Length: $5";

# "Maximum Internal Length: $6";

# "Extend Internal: $7";

# "Length of Start/End Oligomers: $8";

# "Mismatches of Start/End Oligomers: $9";

# "Minimum LTR Length: ${10}";

# "Maximum LTR Length: ${11}";

# "Species code: ${12}";

# "LTR retrotransposon lineage: ${13}"

# "BLAST FULL-LENGTH Coverage (Minimum value): ${14}";

# "BLAST Minimum Signature Coordinate: ${15}"

# "BLAST Maximum Signature Coordinate: ${16}"

# "BLAST PBS Signature Mismatches: ${17}"

# "BLAST PPT Signature Mismatches: ${18}"

# "BLAST Mismatches of Start/End Oligomers: ${19}";

# "BLAST SOLO-LTR Coverage (Minimum value): ${20}";

# "ORF HMM DATABASE: ${21}";

# "HMMSCAN E-VALUE: ${22}";


# Do not forget to unzip the genome fasta file before executing the below command

# Comment out one of the following running commands of the ATHILAfinder based on your operating system


# running command for Linux

bash ATHILAfinder.sh $2 ATHILA_PBSjunction.txt ATHILA_PPTjunction.txt 11000 2000 11000 2000 20 5 1000 2000 Atha Athila 0.80 500 2500 5 5 19 0.90 orfis.Ty3.updated.hmmdb 0.01


# running command for Mac

# script LOGFILE_Atha.txt bash ATHILAfinder.sh 100620.Chr_scaffolds.fa 11000 2000 11000 2000 20 4 1000 2000 Atha 0.90 500 2500 5 5 19 0.98 orfis.Ty3.updated.hmmdb 0.01


exit

