#!/usr/bin/env python
from __future__ import print_function
import sys
import os

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        new_header_count = 1
        current_header = None

        for line in infile:
            line = line.strip()

            if line.startswith('>'):
                # Header line
                current_header = '>Chr{0}'.format(new_header_count)
                new_header_count += 1
                outfile.write(current_header + '\n')
            else:
                # Sequence line
                outfile.write(line.upper() + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_fasta(input_file, output_file)

    print("Conversion complete. Output saved to {0}".format(output_file))
