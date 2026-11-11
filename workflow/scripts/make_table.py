#!/usr/bin/env python3
"""Collect compare_gb.py score files into one Markdown table.
Usage: make_table.py out.md results/score/<variant>/<reference>.txt ..."""
import os
import sys

ORDER = {"raw": 0, "norm": 1}


def read(path):
    rows = {}
    for line in open(path):
        f = line.split()
        if len(f) == 8 and f[0] in ("CDS", "tRNA", "rRNA", "ALL"):
            rows[f[0]] = [int(x) for x in f[1:6]] + [float(f[6]), float(f[7])]
    return rows


def cell(rows, t):
    r = rows.get(t)
    return "n/a" if not r else "%d/%d" % (r[2], r[0])


entries = []
for path in sys.argv[2:]:
    variant = os.path.basename(os.path.dirname(path))
    ref = os.path.splitext(os.path.basename(path))[0]
    entries.append((ref, variant, read(path)))
entries.sort(key=lambda e: (e[0], ORDER.get(e[1], 9)))

lines = ["| reference | names | CDS found | tRNA found | rRNA found | all: sensitivity | all: precision | all: exact start |",
         "|---|---|---|---|---|---|---|---|"]
for ref, variant, rows in entries:
    a = rows.get("ALL")
    sens = "n/a" if not a else "%.1f%%" % (100 * a[5])
    prec = "n/a" if not a else "%.1f%%" % (100 * a[6])
    exact = "n/a" if not a else "%d/%d" % (a[3], a[0])
    lines.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
        ref.replace("_", " "), "as downloaded" if variant == "raw" else "normalized",
        cell(rows, "CDS"), cell(rows, "tRNA"), cell(rows, "rRNA"), sens, prec, exact))
open(sys.argv[1], "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
