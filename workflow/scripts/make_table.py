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
        if len(f) in (8, 11) and f[0] in ("CDS", "tRNA", "rRNA", "ALL"):
            r = [int(x) for x in f[1:6]] + [float(f[6]), float(f[7])]
            r += [int(f[8]), float(f[9]), float(f[10])] if len(f) == 11 else [0, float(f[6]), float(f[7])]
            rows[f[0]] = r
    return rows


def cell(rows, t, near=False):
    r = rows.get(t)
    return "n/a" if not r else "%d/%d" % (r[2] + (r[7] if near else 0), r[0])


entries = []
for path in sys.argv[2:]:
    variant = os.path.basename(os.path.dirname(path))
    ref = os.path.splitext(os.path.basename(path))[0]
    entries.append((ref, variant, read(path)))
entries.sort(key=lambda e: (e[0], ORDER.get(e[1], 9)))

lines = ["| reference | names | CDS found | tRNA found | tRNA found (+-5 nt) | rRNA found | all: sensitivity | all: sensitivity (+-5 nt) | all: precision | all: exact start |",
         "|---|---|---|---|---|---|---|---|---|---|"]
for ref, variant, rows in entries:
    a = rows.get("ALL")
    sens = "n/a" if not a else "%.1f%%" % (100 * a[5])
    sens5 = "n/a" if not a else "%.1f%%" % (100 * a[8])
    prec = "n/a" if not a else "%.1f%%" % (100 * a[6])
    exact = "n/a" if not a else "%d/%d" % (a[3], a[0])
    lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        ref.replace("_", " "), "as downloaded" if variant == "raw" else "normalized",
        cell(rows, "CDS"), cell(rows, "tRNA"), cell(rows, "tRNA", near=True), cell(rows, "rRNA"),
        sens, sens5, prec, exact))
open(sys.argv[1], "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
