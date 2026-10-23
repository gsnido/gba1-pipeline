#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 <input_id_list.txt> <fastq_dir>" >&2
    exit 1
}

if [[ $# -lt 2 || $# -gt 2 ]]; then
    usage
fi

INPUT_FILE="$1"
FASTQ_DIR="$2"

[[ -f "$INPUT_FILE" ]]      || { echo "Error: input file '$INPUT_FILE' not found." >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]]       || { echo "Error: fastq directory '$FASTQ_DIR' not found." >&2; exit 1; }

sanitize_id() {
    local id="$1"
    id="${id//-/_}"
    id="$(echo "$id" | tr -s ' \t' '_')"
    echo "$id"
}

while IFS= read -r raw_id || [[ -n "$raw_id" ]]; do
    [[ -z "$raw_id" ]] && continue
    sanitized_id="$(sanitize_id "$raw_id")"

    matches=()
    while IFS= read -r -d '' f; do
        matches+=("$f")
    done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "${sanitized_id}_*" -print0)

    n_matches=${#matches[@]}

    if [[ "$n_matches" -ne 2 ]]; then
        echo "WARNING: expected exactly 2 matches for ID '${raw_id}' (sanitized: '${sanitized_id}'), found ${n_matches}. Skipping." >&2
        printf '%s\n' "${matches[@]}" >&2
        continue
    fi

    r1=""
    r2=""
    for f in "${matches[@]}"; do
        base="$(basename "$f")"
        if [[ "$base" == *_R1_* || "$base" == *_R1.* ]]; then
            r1="$f"
        elif [[ "$base" == *_R2_* || "$base" == *_R2.* ]]; then
            r2="$f"
        fi
    done

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: could not unambiguously assign R1/R2 for ID '${raw_id}' (sanitized: '${sanitized_id}'). Matches were:" >&2
        printf '%s\n' "${matches[@]}" >&2
        continue
    fi

    #echo "ID '${raw_id}' found!"
    #echo " - R1: $r1"
    #echo " - R2: $r2"

done < "$INPUT_FILE"
