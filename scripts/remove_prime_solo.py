import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
remove_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file
# new_out_file = result_file + "_sl.fasta"

remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(result_file, "w") as f:
    for seq in fasta_sequences:
        nam = seq.id
        nuc = str(seq.seq)
        if nam not in remove and len(nuc) > 0:
            SeqIO.write([seq], f, "fasta")


# try:    
#     with open(new_out_file, "w") as newFasta:
#         for fasta in fasta_sequences:
#             print(fasta)
#             try:
#                 header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
#             except ValueError:
#                 print(fasta)
#             header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
#             sequence = sequence.replace("\n","") + "\n" # Replace newlines in sequence, remember to add one to the end.
#             newFasta.write(header + sequence)
#         newFasta.close()
# except IOError:
#     print("Failed to open " + inFile)
#     exit(1)