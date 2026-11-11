#!/usr/bin/env bash
# Usage: run_pga.sh PGA.pl REF_DIR TARGET.fasta OUT.gb LOG
# Runs PGA in a private temp directory and accepts the result only if it contains
# at least MIN_CDS annotated CDS features (default 30); otherwise it retries
# (PGA_TRIES, default 5) and finally fails. An unannotated file must never be scored.
set -u
pga=$(readlink -f "$1"); refdir=$(readlink -f "$2"); fasta=$(readlink -f "$3")
mkdir -p "$(dirname "$4")" "$(dirname "$5")"
out="$(readlink -f "$(dirname "$4")")/$(basename "$4")"
log="$(readlink -f "$(dirname "$5")")/$(basename "$5")"
min=${MIN_CDS:-30}; tries=${PGA_TRIES:-5}
: > "$log"
for attempt in $(seq 1 "$tries"); do
    work=$(mktemp -d)
    mkdir -p "$work/ref" "$work/target" "$work/out"
    cp "$refdir"/*.gb "$work/ref/"
    cp "$fasta" "$work/target/"
    (cd "$work/out" && perl "$pga" -r ../ref -t ../target > run.log 2>&1)
    gbf="$work/out/gb/$(basename "$fasta" .fasta).gb"
    n=$(grep -c '^     CDS' "$gbf" 2>/dev/null); n=${n:-0}
    echo "attempt $attempt: $n CDS features" >> "$log"
    if [ "$n" -ge "$min" ]; then
        cp "$gbf" "$out"
        cp "$work/out/gb/warning.log" "${out%.gb}.warning.log" 2>/dev/null || true
        rm -rf "$work"
        exit 0
    fi
    tail -5 "$work/out/run.log" >> "$log"
    rm -rf "$work"
done
echo "PGA gave fewer than $min CDS in $tries attempts; see $log" >&2
exit 1
