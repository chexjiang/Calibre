#!/usr/bin/env python3
"""Compute per-cell fragment QC in one streaming pass.

The Satpathy fragment file is very large.  Signac's QC helpers may scan it
several times, so this small helper keeps only the requested barcodes and
computes total fragments, nucleosome signal, and a reproducible TSS-proximity
ratio in one pass without materialising the full fragments table.
"""
import bisect
import gzip
import subprocess
import sys


def main():
    if len(sys.argv) != 5:
        raise SystemExit("usage: fragment_qc.py FRAGMENTS SELECTED_BARCODES TSS_TSV OUTPUT")
    fragments, selected_file, tss_file, output = sys.argv[1:]
    cells = {line.strip() for line in open(selected_file) if line.strip()}
    tss = {}
    for line in open(tss_file):
        chrom, pos = line.rstrip("\n").split("\t")[:2]
        tss.setdefault(chrom, []).append(int(pos))
    for chrom in tss:
        tss[chrom].sort()

    # total fragments, nucleosome-sized fragments, nucleosome-free fragments,
    # cut sites in +/-100 bp and +/-1000 bp around annotated TSSs.
    stats = {c: [0, 0, 0, 0, 0] for c in cells}
    # pigz uses all available cores for decompression and is substantially
    # faster than Python's single-threaded gzip reader on this file.
    stream = subprocess.Popen(["pigz", "-dc", fragments], stdout=subprocess.PIPE,
                              text=True, encoding="utf-8", errors="replace")
    with stream.stdout as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, start, end, cell, count = fields[:5]
            if cell not in stats:
                continue
            start = int(start)
            end = int(end)
            count = int(count)
            length = end - start
            s = stats[cell]
            s[0] += count
            if 147 <= length <= 294:
                s[1] += count
            elif length < 147:
                s[2] += count
            positions = tss.get(chrom)
            if not positions:
                continue
            # Count both cut sites.  A fragment can contribute twice, as in
            # the standard TSS insertion-count calculation.
            for cut in (start + 1, end):
                i = bisect.bisect_left(positions, cut)
                nearby = []
                if i < len(positions):
                    nearby.append(positions[i])
                if i:
                    nearby.append(positions[i - 1])
                if any(abs(cut - p) <= 100 for p in nearby):
                    s[3] += count
                if any(abs(cut - p) <= 1000 for p in nearby):
                    s[4] += count
    return_code = stream.wait()
    if return_code != 0:
        raise SystemExit(f"pigz failed with status {return_code}")

    with open(output, "w") as out:
        out.write("barcode\tfragment_count\tnucleosome_fragments\tnfr_fragments\tcut_core\tcut_window\n")
        for cell in sorted(stats):
            out.write(cell + "\t" + "\t".join(str(x) for x in stats[cell]) + "\n")


if __name__ == "__main__":
    main()
