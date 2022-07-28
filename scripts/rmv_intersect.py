import sys

in_file1=sys.argv[1]
in_file2=sys.argv[2]
out_file=sys.argv[3]

lines_to_remove = None
with open(in_file1) as f0:
    lines_to_remove = {line.rstrip() for line in f0.readlines()}

remaining_lines = None
with open(in_file2) as f1:
    remaining_lines = {line.rstrip() for line in f1.readlines()} - lines_to_remove

with open(out_file, "w") as f2:
    for line in remaining_lines:
        f2.write(line + "\n")