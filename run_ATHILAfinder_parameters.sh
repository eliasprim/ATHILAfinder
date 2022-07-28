# Unzip the zipped assembly file

gzip -d Col-CEN_v1.2.fasta.gz

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

# "BLAST SOLO-LTR Coverage (Minimum value): ${16}";

# "ORF HMM DATABASE: ${17}";


# Do not forget to unzip the genome fasta file before execute the below command

script -c "bash ATHILAfinder.sh Col-CEN_v1.2.fasta 11000 2000 11000 2000 20 4 1000 2000 Atha 0.90 500 2500 5 5 0.98 orfis.Ty3.updated.hmmdb" LOGFILE_Atha.txt

exit