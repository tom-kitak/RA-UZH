#!/usr/bin/env bash
set -euo pipefail

# --- Settings ---
RELEASE=115
SPECIES=homo_sapiens
ASSEMBLY=GRCh38
OUTDIR=${1:-ensembl_release_${RELEASE}}

# --- URLs ---
BASE="https://ftp.ensembl.org/pub/release-${RELEASE}"
GTF_URL="${BASE}/gtf/${SPECIES}/Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf.gz"
DNA_BASE="${BASE}/fasta/${SPECIES}/dna"

# --- Prep output dir ---
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# --- Show the GTF URL ---
echo "GTF file URL:"
echo "${GTF_URL}"
echo

# --- Download GTF ---
echo "Downloading GTF file..."
wget -q --show-progress -c -N "${GTF_URL}"

# --- Download chromosomes (exclude MT) ---
echo
echo "Downloading chromosomes (1–22, X, Y)..."
chroms=($(seq 1 22) X Y)
for c in "${chroms[@]}"; do
  url="${DNA_BASE}/Homo_sapiens.${ASSEMBLY}.dna.chromosome.${c}.fa.gz"
  echo "  -> ${url}"
  wget -q --show-progress -c -N "${url}"
done

# --- Extract all .gz files ---
echo
echo "Extracting all .gz files..."
for f in *.gz; do
  echo "  -> Decompressing ${f}"
  gunzip -f "${f}"   # -f = force overwrite if output exists
done

# --- Remove MT data from GTF in place ---
echo
echo "Removing MT annotations from GTF..."
GTF_FILE="Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf"
if grep -q -w '^MT' "${GTF_FILE}"; then
  tmpfile=$(mktemp)
  grep -v -w '^MT' "${GTF_FILE}" > "${tmpfile}"
  mv "${tmpfile}" "${GTF_FILE}"
  echo "MT lines removed."
else
  echo "No MT lines found; nothing to remove."
fi

echo
echo "All files downloaded, extracted, and cleaned."
echo "Directory: ${OUTDIR}"
