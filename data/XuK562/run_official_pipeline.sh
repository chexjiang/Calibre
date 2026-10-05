#!/usr/bin/env bash
set -euo pipefail

# STARsolo opens many per-cell temporary files while sorting BAM output.
# The default login limit (often 1024) is too low for the K562 run.
ulimit -n 65535 2>/dev/null || true

root="$(cd "$(dirname "$0")/../.." && pwd)"
data_dir="$root/data/XuK562"
sra_dir="$data_dir/sra"
fastq_dir="$data_dir/official_pipeline/fastq"
out="$data_dir/official_pipeline"
mkdir -p "$fastq_dir" "$out/rna" "$out/atac"
mkdir -p "$data_dir/tmp"

sra_tool="${SRA_TOOL:-/tmp/sratoolkit.3.4.1-ubuntu64/bin}"
star="${STAR_BIN:-/tmp/STAR-2.7.9a/bin/Linux_x86_64_static/STAR}"
chromap="${CHROMAP_BIN:-/tmp/chromap/chromap}"
star_index="${STAR_INDEX:-$root/reference/XuK562/star_index}"
chromap_index="${CHROMAP_INDEX:-$root/reference/XuK562/chromap_index/hg38.index}"
whitelist_dir="$root/reference/XuK562/10x_multiome"

for f in "$sra_tool/fasterq-dump" "$star" "$chromap" "$star_index/Genome" "$chromap_index" \
         "$whitelist_dir/gex_737K-arc-v1.txt" "$whitelist_dir/atac_737K-arc-v1_rc.txt"; do
  [[ -e "$f" ]] || { echo "Missing required input: $f" >&2; exit 1; }
done

if [[ ! -s "$fastq_dir/ERR9847049_1.fastq" ||
      ! -s "$fastq_dir/ERR9847049_2.fastq" ||
      ! -s "$fastq_dir/ERR9847049_3.fastq" ]]; then
  rm -f "$fastq_dir"/ERR9847049_*.fastq
  "$sra_tool/fasterq-dump" --split-files --include-technical --threads 32 --temp "$data_dir/tmp" \
    --outdir "$fastq_dir" "$sra_dir/ERR9847049/ERR9847049.sra"
fi
if [[ ! -s "$fastq_dir/ERR9847050_1.fastq" ||
      ! -s "$fastq_dir/ERR9847050_2.fastq" ]]; then
  rm -f "$fastq_dir"/ERR9847050_*.fastq
  "$sra_tool/fasterq-dump" --split-files --threads 32 --temp "$data_dir/tmp" \
    --outdir "$fastq_dir" "$sra_dir/ERR9847050/ERR9847050.sra"
fi

mkdir -p "$out/rna/outs" "$out/atac/outs"
if [[ ! -s "$out/rna/outs/Solo.out/Gene/filtered/matrix.mtx" ]]; then
  "$star" --genomeDir "$star_index" \
    --outFileNamePrefix "$out/rna/outs/" --readFilesIn \
    "$fastq_dir/ERR9847050_2.fastq" "$fastq_dir/ERR9847050_1.fastq" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloType CB_UMI_Simple \
    --clip3pNbases 115 --soloCBwhitelist "$whitelist_dir/gex_737K-arc-v1.txt" \
    --soloCellFilter EmptyDrops_CR --runThreadN 64 \
    --outSAMtype None
fi

if [[ ! -s "$out/atac/outs/fragments.atac.tsv.gz" || ! -s "$out/atac/outs/fragments.atac.tsv.gz.tbi" ]]; then
  if [[ ! -s "$out/atac/outs/fragments.atac.tsv" ]]; then
    "$chromap" --preset atac -t 64 -x "$chromap_index" \
      -r "$root/reference/XuK562/hg38.fa" \
      -1 "$fastq_dir/ERR9847049_2.fastq" -2 "$fastq_dir/ERR9847049_3.fastq" \
      -b "$fastq_dir/ERR9847049_1.fastq" \
      --barcode-whitelist "$whitelist_dir/atac_737K-arc-v1_rc.txt" \
      -o "$out/atac/outs/fragments.atac.tsv"
  fi
  bgzip -@ 64 -f "$out/atac/outs/fragments.atac.tsv"
  tabix -p bed "$out/atac/outs/fragments.atac.tsv.gz"
fi

# 10x Multiome uses distinct GEX and ATAC barcode sequences.  They are paired
# by line number in the two ARC v1 whitelists; translate ATAC barcodes to GEX
# barcodes before constructing the joint object.
if [[ ! -s "$out/atac/outs/fragments.tsv.gz" || ! -s "$out/atac/outs/fragments.tsv.gz.tbi" ]]; then
  translation="$out/atac/arc_barcode_translation.tsv"
  if [[ ! -s "$translation" ]]; then
    paste "$whitelist_dir/atac_737K-arc-v1_rc.txt" \
          "$whitelist_dir/gex_737K-arc-v1.txt" > "$translation"
  fi
  zcat "$out/atac/outs/fragments.atac.tsv.gz" \
    | awk -v mapfile="$translation" '
        BEGIN { while ((getline < mapfile) > 0) { map[$1] = $2 } close(mapfile) }
        { if ($4 in map) $4 = map[$4]; print }
      ' OFS="\t" \
    | bgzip -@ 64 > "$out/atac/outs/fragments.tsv.gz"
  tabix -p bed "$out/atac/outs/fragments.tsv.gz"
fi

# Call aggregate peaks from an even genomic subsample of fragments.  MACS2's
# paired-end parser requires a seekable input, so this avoids materializing
# hundreds of millions of pseudo-reads while retaining genome-wide coverage.
peak_dir="$out/atac/macs2_full"
mkdir -p "$peak_dir"
if [[ ! -s "$peak_dir/aggregate_peaks.narrowPeak" ]]; then
  sample_bedpe="$out/atac/subset_fragments.bedpe"
  zcat "$out/atac/outs/fragments.atac.tsv.gz" \
    | awk 'BEGIN{OFS="\t"} NR % 60 == 0 {print $1,$2,$3,".",1,"+"}' \
    > "$sample_bedpe"
  macs2 callpeak -t "$sample_bedpe" -g hs -f BEDPE \
    --nomodel --keep-dup all -B --SPMR --outdir "$peak_dir" -n aggregate
fi
cut -f 1-4 "$peak_dir/aggregate_peaks.narrowPeak" \
  | sed '/chrM/d' \
  | bedtools intersect -a - -b "$root/reference/XuK562/hg38-blacklist.bed.gz" -v \
  | sort -k1,1 -k2,2n > "$peak_dir/aggregate_peaks_formatted.bed"

echo "RNA and ATAC alignment completed. Run assemble_official_xuk562.R next."
