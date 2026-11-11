#!/usr/bin/env python3
"""Rename genes in a RefSeq plastome GenBank file so that PGA recognises them.

PGA (Qu et al. 2019) expects
  * rRNA genes named rrn16 / rrn23 / rrn4.5 / rrn5   (RefSeq: rrn16S, rrn16s, ...)
  * tRNA genes named with an anticodon, e.g. trnQ-UUG (RefSeq often: trnQ)
Genes that do not follow this convention are silently skipped by PGA.

This script
  1. strips the trailing S/s from rrn names;
  2. adds the anticodon to bare tRNA names. The anticodon is read from the
     spliced tRNA sequence (anticodon loop: Y-U-[NNN]-purine) and checked against
     the amino acid in /product (anticodon -> codon -> amino acid).

Usage: normalize_reference.py in.gb out.gb report.tsv   (keep the report OUT of the reference folder)
"""
import re
import sys

BASES = "tcag"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AAS[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "fMet": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T",
       "Trp": "W", "Tyr": "Y", "Val": "V"}
COMP = str.maketrans("acgtn", "tgcan")
RRN = re.compile(r"^rrn(\d+(?:\.\d+)?)[sS]$")
LONG = re.compile(r"^tRNA-([A-Za-z]+) \(([ACGUTacgut]{3})\)$")
FEAT = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
GENE = re.compile(r'/gene="([^"]*)"')


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse_loc(loc):
    strand = "-" if "complement" in loc else "+"
    nums = re.findall(r"(\d+)\.\.(\d+)", loc.replace("<", "").replace(">", ""))
    return strand, sorted((int(a), int(b)) for a, b in nums)


def mature(seq, loc):
    strand, exons = parse_loc(loc)
    if not exons:
        return None
    g = "".join(seq[a - 1:b] for a, b in exons)
    return revcomp(g) if strand == "-" else g


def find_anticodon(m, aa1):
    cands = []
    for s in range(28, min(42, len(m) - 4)):
        if m[s - 1] == "t" and m[s - 2] in "ct" and m[s + 3] in "ag":
            ac = m[s:s + 3]
            aa = CODON.get(revcomp(ac))
            ok = aa == aa1 or (aa1 == "I" and ac == "cat")  # lysidine tRNA-Ile(CAU)
            cands.append((not ok, abs(s - 33), ac))
    if not cands:
        return None, "failed"
    cands.sort()
    bad, _, ac = cands[0]
    return ac, ("unverified" if bad else "ok")


def main(inp, outp, rep):
    lines = open(inp).read().split("\n")
    fi = next(i for i, l in enumerate(lines) if l.startswith("FEATURES"))
    oi = next(i for i, l in enumerate(lines) if l.startswith("ORIGIN"))
    body = "\n".join(lines[oi + 1:]).split("//")[0].lower()
    seq = "".join(re.findall(r"[a-z]+", body))

    blocks = []
    i = fi + 1
    while i < oi:
        m = FEAT.match(lines[i])
        if m and not m.group(1).startswith("/"):
            j = i + 1
            while j < oi and not FEAT.match(lines[j]):
                j += 1
            loc = m.group(2)
            for k in range(i + 1, j):
                if lines[k].strip().startswith("/"):
                    break
                loc += lines[k].strip()
            blocks.append({"key": m.group(1), "a": i, "b": j, "loc": loc})
            i = j
        else:
            i += 1

    def gene_of(b):
        for k in range(b["a"], b["b"]):
            g = GENE.search(lines[k])
            if g:
                return g.group(1)
        return None

    def rename(b, old, new):
        for k in range(b["a"], b["b"]):
            lines[k] = lines[k].replace('/gene="%s"' % old, '/gene="%s"' % new)

    def product_aa(b, name):
        for k in range(b["a"], b["b"]):
            p = re.search(r'/product="tRNA-([A-Za-z]+)', lines[k])
            if p and p.group(1) in AA3:
                return AA3[p.group(1)]
        n = name[3:]
        return "M" if n.startswith("fM") else (n[:1] if n[:1].isalpha() else None)

    report = []
    for n, b in enumerate(blocks):
        name = gene_of(b)
        if not name:
            continue
        m = RRN.match(name)
        if m:
            new = "rrn" + m.group(1)
            rename(b, name, new)
            report.append((b["key"], name, new, b["loc"], "rrn suffix removed", ""))
            continue
        if b["key"] != "tRNA":
            continue
        lg = LONG.match(name)
        if lg and lg.group(1) in AA3:
            letter = "fM" if lg.group(1) == "fMet" else AA3[lg.group(1)]
            new = "trn%s-%s" % (letter, lg.group(2).upper().replace("T", "U"))
            rename(b, name, new)
            if n > 0 and blocks[n - 1]["key"] == "gene" and gene_of(blocks[n - 1]) == name:
                rename(blocks[n - 1], name, new)
            report.append(("tRNA", name, new, b["loc"], "converted from long name", ""))
            continue
        if not name.startswith("trn"):
            report.append(("tRNA", name, name, b["loc"], "WARNING: not a trn name", ""))
            continue
        if "-" in name:
            continue
        mat = mature(seq, b["loc"])
        aa1 = product_aa(b, name)
        if mat is None:
            report.append(("tRNA", name, name, b["loc"], "FAILED: location", ""))
            continue
        ac, status = find_anticodon(mat, aa1)
        if ac is None:
            report.append(("tRNA", name, name, b["loc"], "FAILED: no anticodon loop", mat))
            continue
        new = "%s-%s" % (name, ac.upper().replace("T", "U"))
        rename(b, name, new)
        if n > 0 and blocks[n - 1]["key"] == "gene" and gene_of(blocks[n - 1]) == name:
            rename(blocks[n - 1], name, new)
        report.append(("tRNA", name, new, b["loc"], status,
                       "%s->%s" % (revcomp(ac).upper(), CODON.get(revcomp(ac)))
                       + (" " + mat if status != "ok" else "")))

    open(outp, "w").write("\n".join(lines))
    with open(rep, "w") as r:
        r.write("feature\told\tnew\tlocation\tstatus\tcodon->aa\n")
        for row in report:
            r.write("\t".join(row) + "\n")
    cnt = {}
    for row in report:
        cnt[row[4]] = cnt.get(row[4], 0) + 1
    for row in report:
        print("%-5s %-9s -> %-11s %-12s %s" % (row[0], row[1], row[2], row[4], row[5]))
    print("summary:", cnt)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
