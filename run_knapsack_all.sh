#!/usr/bin/env bash
# Run ./memobound on every file in Instance/knapsack; write full log to resultat.txt
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT" || exit 1

BIN="$ROOT/build/memoscipl"
OUT="$ROOT/Instance/knapsack/resultatL.txt"
INST="$ROOT/Instance/knapsack"

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found or not executable (build first, run from repo root)" >&2
    exit 1
fi
if [[ ! -d "$INST" ]]; then
    echo "error: directory not found: $INST" >&2
    exit 1
fi

: >"$OUT"

while IFS= read -r -d '' f; do
    echo "Instance : $(basename "$f")" >&2
    {
        echo ""
        echo "========== $f =========="
        echo "$BIN $f -l=1"
        "$BIN" "$f" -l=1
        ec=$?
        if [[ $ec -ne 0 ]]; then
            echo "[exit code: $ec]"
        fi
    } >>"$OUT" 2>&1
done < <(find "$INST" -maxdepth 1 -type f -name '*.lp' -print0 | sort -z)

echo "Wrote $OUT"
