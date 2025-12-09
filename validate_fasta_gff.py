#!/usr/bin/env python3
import sys

def load_fasta_lengths(fasta_file):
    lengths = {}
    with open(fasta_file) as f:
        seqid = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if seqid:
                    lengths[seqid] = length
                seqid = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
        if seqid:
            lengths[seqid] = length
    return lengths

def validate_gff_coords(gff_file, fasta_lengths):
    missing_seqids = set()
    invalid_coords = []
    with open(gff_file) as f:
        for lineno, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                print(f"Warning: Line {lineno} does not have 9 columns, skipping")
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if seqid not in fasta_lengths:
                missing_seqids.add(seqid)
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                print(f"Warning: Line {lineno} has non-integer start/end: {start}, {end}")
                continue
            seq_len = fasta_lengths[seqid]
            if start_i < 1 or end_i > seq_len or start_i > end_i:
                invalid_coords.append((lineno, seqid, start_i, end_i, seq_len))
    return missing_seqids, invalid_coords

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 validate_fasta_gff.py <input.gff> <input.fasta>")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]

    print(f"Loading FASTA lengths from {fasta_file} ...")
    fasta_lengths = load_fasta_lengths(fasta_file)
    print(f"Found {len(fasta_lengths)} sequences in FASTA.")

    print(f"Validating GFF coordinates in {gff_file} ...")
    missing_seqids, invalid_coords = validate_gff_coords(gff_file, fasta_lengths)

    if missing_seqids:
        print(f"ERROR: The following seqids are in GFF but missing in FASTA ({len(missing_seqids)}):")
        for s in sorted(missing_seqids):
            print(f"  {s}")
    else:
        print("All seqids in GFF are present in FASTA.")

    if invalid_coords:
        print(f"ERROR: Found {len(invalid_coords)} features with invalid coordinates:")
        for lineno, seqid, start_i, end_i, seq_len in invalid_coords[:20]:
            print(f"  Line {lineno}: seqid={seqid}, start={start_i}, end={end_i}, seq_len={seq_len}")
        if len(invalid_coords) > 20:
            print(f"  ... and {len(invalid_coords)-20} more")
    else:
        print("All feature coordinates are valid within FASTA sequence lengths.")

if __name__ == "__main__":
    main()

