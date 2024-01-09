#!/bin/bash

# Unzip the zipped assembly file

gzip -d Col-CEN_v1.2.fasta.gz

# Pre-process assembly file

python treat_fasta.py Col-CEN_v1.2.fasta


# "Assembly: $1";

# "Searching Regions (bedtools window): $2";

# "Minimum Internal Length: $3";

# "Maximum Internal Length: $4";

# "Extend Internal: $5";

# "Length of Start/End Oligomers: $6";

# "Mismatches of Start/End Oligomers: $7";

# "Minimum LTR Length: $8";

# "Maximum LTR Length: $9";

# "Species code: ${10}";

# "BLAST FULL-LENGTH Coverage (Minimum value): ${11}";

# "BLAST Minimum Signature Coordinate: ${12}"

# "BLAST Maximum Signature Coordinate: ${13}"

# "BLAST PBS Signature Mismatches: ${14}"

# "BLAST PPT Signature Mismatches: ${15}"

# "Mismatches of Start/End Oligomers BLAST: ${16}";

# "BLAST SOLO-LTR Coverage (Minimum value): ${17}";

# "ORF HMM DATABASE: ${18}";

# "HMMSCAN E-VALUE: ${19}";


# Do not forget to unzip the genome fasta file before executing the below command

# Comment out one of the following running commands of the ATHILAfinder based on your operating system


# running command for Linux

script -c "bash ATHILAfinder.sh Col-CEN_v1.2_clean.fasta 11000 2000 11000 2000 20 4 1000 2000 Atha 0.90 500 2500 5 5 15 0.98 orfis.Ty3.updated.hmmdb 0.01" LOGFILE_Atha.txt


# running command for Mac

# script LOGFILE_Atha.txt bash ATHILAfinder.sh Col-CEN_v1.2.fasta 11000 2000 11000 2000 20 4 1000 2000 Atha 0.90 500 2500 5 5 15 0.98 orfis.Ty3.updated.hmmdb 0.01


exit

