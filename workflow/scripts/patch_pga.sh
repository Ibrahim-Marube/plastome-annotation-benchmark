#!/usr/bin/env bash
# Usage: patch_pga.sh PGA.pl PGA.patched.pl
# PGA builds temp-file names from 14 random characters; its character list contains '$'.
# Those names go straight into shell commands, where '$$' becomes a PID and '$name'
# expands to nothing, so makeblastdb/blastn receive a different file name than the one
# PGA wrote. The run then "succeeds" with an unannotated output. We remove '$' from the list.
set -euo pipefail
cp "$1" "$2"
perl -pi -e 's/\x27\$\x27,// if /^my \@a =/' "$2"
if grep -q "^my @a = .*'[$]'" "$2"; then
    echo "patch failed: '\$' still in the character list" >&2
    exit 1
fi
changed=$(diff "$1" "$2" | grep -c '^>' || true)
if [ "$changed" -ne 1 ]; then
    echo "patch changed $changed lines, expected exactly 1" >&2
    exit 1
fi
echo "patched $2 (1 line changed)"
