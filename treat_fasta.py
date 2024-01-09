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
                current_header = f'>Chr{new_header_count}'
                new_header_count += 1
                outfile.write(current_header + '\n')
            else:
                # Sequence line
                outfile.write(line.upper() + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fasta")
        sys.exit(1)

    input_file = sys.argv[1]

    # Remove the extension from the input file
    input_filename_without_extension = os.path.splitext(input_file)[0]

    # Generate the output file name by appending "_clean" before the file extension
    output_file = f"{input_filename_without_extension}_clean.fasta"

    process_fasta(input_file, output_file)

    print(f"Conversion complete. Output saved to {output_file}")