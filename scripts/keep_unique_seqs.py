import sys
from Bio import SeqIO

in_file=sys.argv[1]
out_file=sys.argv[2]

seen = set()
records = []

for record in SeqIO.parse(in_file, "fasta"):  
    if record.seq not in seen:
        seen.add(record.seq)
        records.append(record)


#writing to a fasta file
SeqIO.write(records, out_file, "fasta") 