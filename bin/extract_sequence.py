#!/usr/bin/env python3

import subprocess
import sys

taxid_to_find = sys.argv[1]
output_file = sys.argv[2]
source_fna_path = sys.argv[3]

with open(source_fna_path, "r") as source_file:
    header = ""
    seq = ""
    found = False
    nfound = 0
    for line in source_file:
        line = line.strip()
        if line.startswith(">"):
            if found:
                break
            header = line
            if f"|{taxid_to_find}|" in header:
                nfound += 1
                found = True
        elif found:
            seq += line

if found:
    header_to_write = header.replace("(", "-").replace(")", "")
    with open(output_file, "w") as output:
        output.write(header_to_write + "\n" + seq + "\n")
    print(header_to_write, end="")
    subprocess.run(["bwa", "index", output_file])
