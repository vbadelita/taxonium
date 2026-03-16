#!/usr/bin/env python3
"""Check what percent of HA sequences are present in each other segment."""

import glob
import os

testtree = os.path.join(os.path.dirname(__file__), "testtree")


def read_fasta_names(path):
    names = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                names.add(line[1:].strip())
    return names


ha_names = read_fasta_names(os.path.join(testtree, "HA_filtered.fasta"))
print(f"HA: {len(ha_names)} sequences (reference)\n")

segments = sorted(
    os.path.basename(p).replace("_filtered.fasta", "")
    for p in glob.glob(os.path.join(testtree, "*_filtered.fasta"))
    if "HA" not in os.path.basename(p))

for seg in segments:
    path = os.path.join(testtree, f"{seg}_filtered.fasta")
    names = read_fasta_names(path)
    overlap = names & ha_names
    pct = 100 * len(overlap) / len(ha_names)
    print(
        f"{seg}: {len(names)} sequences, {len(overlap)}/{len(ha_names)} in HA ({pct:.1f}%)"
    )
