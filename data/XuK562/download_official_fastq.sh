#!/usr/bin/env bash
set -euo pipefail

out_dir="$(cd "$(dirname "$0")" && pwd)/raw_fastq"
mkdir -p "$out_dir"
chunk_bytes=$((256 * 1024 * 1024))
manifest="$out_dir/.chunk_manifest.tsv"
: > "$manifest"

add_file() {
  local name="$1" url="$2" size="$3"
  local n=$(( (size + chunk_bytes - 1) / chunk_bytes ))
  rm -f "$out_dir/$name"
  for ((i=0; i<n; i++)); do
    local start=$((i * chunk_bytes))
    local end=$((start + chunk_bytes - 1))
    if (( end >= size )); then end=$((size - 1)); fi
    local part
    part=$(printf '%s.part.%05d' "$name" "$i")
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$start" "$end" "$part" >> "$manifest"
  done
}

add_file K562_10x_Multiome_ATAC_S1_L001_I2_001.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847049/K562_10x_Multiome_ATAC_S1_L001_I2_001.fastq.gz 6105125692
add_file K562_10x_Multiome_ATAC_S1_L001_R1_001.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847049/K562_10x_Multiome_ATAC_S1_L001_R1_001.fastq.gz 38644167854
add_file K562_10x_Multiome_ATAC_S1_L001_R2_001.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847049/K562_10x_Multiome_ATAC_S1_L001_R2_001.fastq.gz 38834220493
add_file K562_10x_Multiome_RNA_S1_L001_R1_001.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847050/K562_10x_Multiome_RNA_S1_L001_R1_001.fastq.gz 35880695685
add_file K562_10x_Multiome_RNA_S1_L001_R2_001.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847050/K562_10x_Multiome_RNA_S1_L001_R2_001.fastq.gz 35534657700

download_chunk() {
  local line="$1"
  IFS=$'\t' read -r name url start end part <<< "$line"
  local target="$out_dir/$part"
  local expected=$((end - start + 1))
  if [[ -s "$target" ]] && [[ "$(stat -c%s "$target")" -eq "$expected" ]]; then return 0; fi
  rm -f "$target"
  curl -L --fail --retry 8 --retry-delay 5 --connect-timeout 30 -sS -r "$start-$end" -o "$target" "$url"
  [[ "$(stat -c%s "$target")" -eq "$expected" ]]
}
export out_dir
export -f download_chunk

xargs -P 8 -I{} bash -c 'download_chunk "$1"' _ "{}" < "$manifest"

for name in \
  K562_10x_Multiome_ATAC_S1_L001_I2_001.fastq.gz \
  K562_10x_Multiome_ATAC_S1_L001_R1_001.fastq.gz \
  K562_10x_Multiome_ATAC_S1_L001_R2_001.fastq.gz \
  K562_10x_Multiome_RNA_S1_L001_R1_001.fastq.gz \
  K562_10x_Multiome_RNA_S1_L001_R2_001.fastq.gz; do
  cat "$out_dir/$name.part."* > "$out_dir/$name"
done

rm -f "$manifest" "$out_dir"/*.part.*
echo "Official Xu K562 FASTQ download complete: $out_dir"
