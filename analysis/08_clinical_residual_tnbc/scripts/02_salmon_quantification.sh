#!/usr/bin/env bash
# Documented project command specification, not an original execution log.
set -euo pipefail
mode=${1:?Use index or quantify}; shift
if ! salmon --version | grep -q '1.10.3'; then printf '%s\n' 'Salmon 1.10.3 is required' >&2; exit 1; fi
if [[ "$mode" == index ]]; then
  transcripts=${1:?GENCODE v36 transcript FASTA}; genome=${2:?GRCh38 full genome FASTA}; out=${3:?Index output directory}
  mkdir -p "$out"
  # Inputs are uncompressed complete FASTA files. Decoys are all genome contigs.
  awk '/^>/{sub(/^>/,""); split($0,a," "); print a[1]}' "$genome" > "$out/decoys.txt"
  cat "$transcripts" "$genome" > "$out/gentrome.fa"
  salmon index -t "$out/gentrome.fa" -d "$out/decoys.txt" -i "$out/index" -k 31 --gencode
elif [[ "$mode" == quantify ]]; then
  index=${1:?Index directory}; read1=${2:?Paired read 1}; read2=${3:?Paired read 2}; out=${4:?Sample output directory}; threads=${5:-8}
  salmon quant -i "$index" -l A -1 "$read1" -2 "$read2" --seqBias --gcBias -p "$threads" -o "$out"
else
  printf '%s\n' 'Use index or quantify' >&2;exit 1
fi
