#!/usr/bin/env python3
"""Score a plastome annotation (e.g. PGA output) against a curated GenBank record.

Features compared: CDS, tRNA, rRNA (pseudogenes in the truth set are skipped).
A feature is identified by (type, strand, 3' end). For each feature type we report
  found      : same strand and same 3' end (stop codon / last base)
  exact      : found, and the 5' end (start) also agrees
  structure  : exact, and every exon boundary agrees (so introns are right too)

Usage: compare_gb.py truth.gb predicted.gb [out_prefix]
"""
import re
import sys
from collections import defaultdict

TYPES = ("CDS", "tRNA", "rRNA")


def split_top(s):
    """Split on commas that are not inside parentheses."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def parse_location(s):
    """Return exons [(start, end, strand)] in biological (5'->3') order."""
    s = s.strip()
    if s.startswith("complement("):
        inner = parse_location(s[len("complement("):-1])
        flip = {"+": "-", "-": "+"}
        return [(a, b, flip[st]) for a, b, st in reversed(inner)]
    if s.startswith(("join(", "order(")):
        body = s[s.index("(") + 1:-1]
        out = []
        for part in split_top(body):
            out.extend(parse_location(part))
        return out
    m = re.match(r"[<>]?(\d+)(?:\.\.[<>]?(\d+))?$", s)
    if not m:
        raise ValueError(f"cannot parse location: {s}")
    a = int(m.group(1))
    b = int(m.group(2) or a)
    return [(a, b, "+")]


def read_features(path):
    feats, in_feat = [], False
    cur = None
    for line in open(path):
        if line.startswith("FEATURES"):
            in_feat = True
            continue
        if line.startswith("ORIGIN"):
            break
        if not in_feat:
            continue
        if line[:5] == "     " and line[5] != " ":
            if cur:
                feats.append(cur)
            key = line[5:21].strip()
            cur = {"type": key, "loc": line[21:].strip(), "quals": {}, "in_loc": True}
        elif cur is not None:
            text = line[21:].strip()
            if text.startswith("/"):
                cur["in_loc"] = False
                k, _, v = text[1:].partition("=")
                cur["quals"][k] = v.strip('"')
            elif cur["in_loc"]:
                cur["loc"] += text
    if cur:
        feats.append(cur)
    return feats


def to_genes(path, truth=False):
    genes = {}
    for f in read_features(path):
        if f["type"] not in TYPES:
            continue
        if truth and ("pseudo" in f["quals"] or "pseudogene" in f["quals"]):
            continue
        exons = parse_location(f["loc"])
        strand = exons[0][2]
        first, last = exons[0], exons[-1]
        start = first[0] if strand == "+" else first[1]
        stop = last[1] if strand == "+" else last[0]
        name = f["quals"].get("gene", f["quals"].get("product", "?"))
        genes[(f["type"], strand, stop)] = {
            "start": start, "stop": stop, "name": name,
            "exons": tuple(sorted((a, b) for a, b, _ in exons)),
        }
    return genes


def main():
    truth_path, pred_path = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else "comparison"
    truth, pred = to_genes(truth_path, truth=True), to_genes(pred_path)
    if len(pred) < 10 or len(truth) < 10:
        sys.exit(f"ERROR: {len(pred)} predicted / {len(truth)} truth features; "
                 f"refusing to score an unannotated file ({pred_path})")
    print(f"{'type':5s} {'truth':>6s} {'pred':>6s} {'found':>6s} {'exact':>6s} {'struct':>6s}   "
          f"sens(found)  prec(found)")
    rows = []
    for t in TYPES + ("ALL",):
        tk = {k for k in truth if t == "ALL" or k[0] == t}
        pk = {k for k in pred if t == "ALL" or k[0] == t}
        found = tk & pk
        exact = {k for k in found if truth[k]["start"] == pred[k]["start"]}
        struct = {k for k in exact if truth[k]["exons"] == pred[k]["exons"]}
        sens = len(found) / len(tk) if tk else float("nan")
        prec = len(found) / len(pk) if pk else float("nan")
        print(f"{t:5s} {len(tk):6d} {len(pk):6d} {len(found):6d} {len(exact):6d} {len(struct):6d}   "
              f"{sens:10.3f}  {prec:10.3f}")
        rows.append((t, len(tk), len(pk), len(found), len(exact), len(struct)))
    for name, keys, src in (("missed", truth.keys() - pred.keys(), truth),
                            ("extra", pred.keys() - truth.keys(), pred)):
        with open(f"{prefix}.{name}.tsv", "w") as out:
            out.write("type\tstrand\tname\tstart\tstop\n")
            for k in sorted(keys, key=lambda k: (k[0], src[k]["name"])):
                g = src[k]
                out.write(f"{k[0]}\t{k[1]}\t{g['name']}\t{g['start']}\t{g['stop']}\n")
    with open(f"{prefix}.start_or_intron_diff.tsv", "w") as out:
        out.write("type\tstrand\tname\ttruth_start\tpred_start\texons_equal\n")
        for k in sorted(truth.keys() & pred.keys()):
            a, b = truth[k], pred[k]
            if a["start"] != b["start"] or a["exons"] != b["exons"]:
                out.write(f"{k[0]}\t{k[1]}\t{a['name']}\t{a['start']}\t{b['start']}\t"
                          f"{a['exons'] == b['exons']}\n")


if __name__ == "__main__":
    main()
