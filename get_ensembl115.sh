#!/usr/bin/env bash
set -euo pipefail

# basic settings
RELEASE=115
SPECIES=homo_sapiens
ASSEMBLY=GRCh38
OUTDIR=${1:-ensembl_release_${RELEASE}}

BASE="https://ftp.ensembl.org/pub/release-${RELEASE}"
GTF_URL="${BASE}/gtf/${SPECIES}/Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf.gz"
DNA_BASE="${BASE}/fasta/${SPECIES}/dna"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "GTF URL: ${GTF_URL}"
echo "Downloading GTF..."
echo "${GTF_URL}"
wget -q --show-progress -c -N "${GTF_URL}"

echo "Downloading chromosomes..."
chroms=($(seq 1 22) X Y)
for c in "${chroms[@]}"; do
  url="${DNA_BASE}/Homo_sapiens.${ASSEMBLY}.dna.chromosome.${c}.fa.gz"
  echo "${url}"
  wget -q --show-progress -c -N "${url}"
done

echo "Extracting .gz files..."
for f in *.gz; do
  gunzip -f "${f}"
done

echo "Removing MT entries from GTF..."
GTF_FILE="Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf"
if grep -q -w '^MT' "${GTF_FILE}"; then
  tmpfile=$(mktemp)
  grep -v -w '^MT' "${GTF_FILE}" > "${tmpfile}"
  mv "${tmpfile}" "${GTF_FILE}"
fi

echo "Renaming FASTA files..."
shopt -s nullglob
for f in Homo_sapiens.${ASSEMBLY}.dna.chromosome.*.fa; do
  newname=$(echo "$f" | sed -E 's/^Homo_sapiens\.GRCh38\.dna\.chromosome\.([^.]+)\.fa$/\1.fa/')
  if [[ "$f" != "$newname" ]]; then
    mv -f "$f" "$newname"
  fi
done
shopt -u nullglob

echo "Done. Files in ${OUTDIR}"
