#!/usr/bin/env python3
"""Extract the sequence of a GenBank record as FASTA (no annotation).
Usage: gb_to_fasta.py in.gb out.fasta   (record name = output file stem)"""
import os
import re
import sys

text = open(sys.argv[1]).read()
if "ORIGIN" not in text:
    sys.exit("no ORIGIN section in " + sys.argv[1])
seq = "".join(re.findall(r"[a-z]+", text.split("ORIGIN")[1].split("//")[0].lower())).upper()
name = os.path.splitext(os.path.basename(sys.argv[2]))[0]
with open(sys.argv[2], "w") as out:
    out.write(">%s\n" % name)
    for i in range(0, len(seq), 70):
        out.write(seq[i:i + 70] + "\n")
print("%s: %d bp" % (name, len(seq)))
