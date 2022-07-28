import sys
import regex as re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

in_file=sys.argv[1]
out_file=in_file + "_ml.fasta"
new_out_file=in_file + "_sl.fasta"


# print(">> Returns the reverse complement of the input fasta")

def make_rc_record(record):
    """Returns a new SeqRecord with the reverse complement sequence."""
    return SeqRecord(seq = record.seq.reverse_complement(), \
                 id = "rc_" + record.id + "_P_element" , \
                 description = "")

records = map(make_rc_record, SeqIO.parse(in_file, "fasta"))
SeqIO.write(records, out_file, "fasta")


# print(">> Opening FASTA file...")
# Reads sequence file list and stores it as a string object. Safely closes file:
try:    
    with open(out_file, "r") as newFile:
        sequences = newFile.read()
        sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
        del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
                         # Del removes this empty element.
        newFile.close()
except IOError:
    print("Failed to open " + inFile)
    exit(1)


# print(">> Converting FASTA file from multiline to single line and writing to file.")
# Conversts multiline fasta to single line. Writes new fasta to file.
try:    
    with open(new_out_file, "w") as newFasta:
        for fasta in sequences:
            try:
                header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
            except ValueError:
                print(fasta)
            header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
            sequence = sequence.replace("\n","") + "\n" # Replace newlines in sequence, remember to add one to the end.
            newFasta.write(header + sequence)
        newFasta.close()
except IOError:
    print("Failed to open " + inFile)
    exit(1)

# print(">> Done!")




# OLD
# output is in multiple line format, run the following command to convert it into single line format

# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  t2t-col.20201227.fasta.F2B.P_fulllength_ml.fasta > t2t-col.20201227.fasta.F2B.P_fulllength_sl.fasta