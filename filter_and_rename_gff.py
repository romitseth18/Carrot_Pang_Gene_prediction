#!/usr/bin/env python3
import sys, re
from collections import defaultdict, OrderedDict

# ------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------

def parse_chr_number(seqid: str):
    m = re.search(r'Chr(\d+)', seqid, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None

def parse_attributes(attr_str):
    attrs = {}
    for kv in attr_str.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            k, v = k.strip(), v.strip()
            if k and v:
                attrs[k] = v
        elif kv.strip():
            print(f"Warning: malformed attribute '{kv}'", file=sys.stderr)
    return attrs

def format_attributes(attrs_dict):
    return ";".join(f"{k}={v}" for k, v in attrs_dict.items())

def get_fasta_lengths(fasta_file):
    lengths = {}
    with open(fasta_file) as f:
        seqid, length = None, 0
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

def get_masked_regions(fasta_file):
    masked = defaultdict(list)
    chrom, start, pos = None, None, 0
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                pos, start = 0, None
                continue
            for base in line:
                if base.upper() == "N":
                    if start is None:
                        start = pos
                else:
                    if start is not None:
                        masked[chrom].append((start, pos))
                        start = None
                pos += 1
        if start is not None:
            masked[chrom].append((start, pos))
    return masked

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(end, last_end))
        else:
            merged.append((start, end))
    return merged

def overlaps(chrom, start, end, masked_dict):
    if chrom not in masked_dict:
        return False
    for s, e in masked_dict[chrom]:
        if end > s and start < e:
            return True
    return False

def parse_fasta_order(fasta_file):
    order = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                order.append(line[1:].split()[0])
    return order

# ------------------------------------------------------------
# Main Logic
# ------------------------------------------------------------

def rename_and_filter_gff_sorted(gff_file, prefix, fasta_file, out_gff, idmap_file):
    fasta_lengths = get_fasta_lengths(fasta_file)
    masked_regions = get_masked_regions(fasta_file)
    for c in masked_regions:
        masked_regions[c] = merge_intervals(masked_regions[c])
    seq_order = parse_fasta_order(fasta_file)
    seq_index = {s: i for i, s in enumerate(seq_order)}

    tag_map = {
        "gene": "gene", "mRNA": "mRNA", "transcript": "trans",
        "CDS": "cds", "exon": "exon", "start_codon": "start",
        "stop_codon": "stop", "tRNA": "tRNA", "intron": "intron"
    }

    gff_by_seq = defaultdict(list)
    headers = []

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                headers.append(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) == 9:
                gff_by_seq[parts[0]].append(line.strip())

    sorted_seqids = sorted(gff_by_seq.keys(), key=lambda x: seq_index.get(x, 1e9))
    contig_map, contig_start = OrderedDict(), 10
    global_gene_counter = 0
    old2new, idmap, features_out = {}, [], []

    for seqid in sorted_seqids:
        def start_pos(l): return int(l.split("\t")[3])
        for line in sorted(gff_by_seq[seqid], key=start_pos):
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            seqid_line, source, ftype, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            if seqid_line not in fasta_lengths:
                continue
            if start > end or end > fasta_lengths[seqid_line]:
                continue

            attrs_dict = parse_attributes(attrs)
            old_id = attrs_dict.get("ID")
            parent_id = attrs_dict.get("Parent")

            chr_num = parse_chr_number(seqid_line)
            if chr_num and 1 <= chr_num <= 9:
                prefix_code = f"{chr_num:02d}"
            else:
                if seqid_line not in contig_map:
                    contig_map[seqid_line] = contig_start
                    contig_start += 1
                prefix_code = f"{contig_map[seqid_line]:02d}"

            # --- Assign hierarchical IDs ---
            if ftype == "gene":
                if overlaps(seqid_line, start, end, masked_regions):
                    continue
                global_gene_counter += 1
                new_id = f"{prefix}_g{prefix_code}{global_gene_counter:08d}"
                old2new[old_id] = new_id
            else:
                tag = tag_map.get(ftype)
                if not tag:
                    continue
                if parent_id and parent_id in old2new:
                    parent_new = old2new[parent_id]
                    new_id = f"{parent_new}_{tag}"
                else:
                    new_id = f"{prefix}_{tag}{prefix_code}{global_gene_counter:08d}"
                old2new[old_id] = new_id

            attrs_dict["ID"] = new_id
            if parent_id:
                if parent_id in old2new:
                    attrs_dict["Parent"] = old2new[parent_id]
                else:
                    attrs_dict["Parent"] = parent_id

            if ftype == "CDS":
                phase = parts[7].strip()
                if phase not in ['0', '1', '2']:
                    phase = '0'
                parts[7] = phase

            parts[8] = format_attributes(attrs_dict)
            parts[1] = source
            features_out.append((ftype, "\t".join(parts)))
            if old_id:
                idmap.append((old_id, new_id))

    # --- Deduplicate gene entries ---
    seen_genes = set()
    final_features = []
    for ftype, line in features_out:
        if ftype == "gene":
            gid = re.search(r'ID=([^;]+)', line)
            if gid:
                gid = gid.group(1)
                if gid in seen_genes:
                    continue
                seen_genes.add(gid)
        final_features.append(line)

    with open(out_gff, "w") as out:
        for h in headers:
            out.write(h)
        for line in final_features:
            out.write(line + "\n")

    with open(idmap_file, "w") as f:
        for old, new in idmap:
            f.write(f"{old}\t{new}\n")

# ------------------------------------------------------------
# Entry Point
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 filter_and_rename_gff_fixed.py <input.gff> <prefix> <fasta> <output.gff>")
        sys.exit(1)

    gff, prefix, fasta, out = sys.argv[1:5]
    rename_and_filter_gff_sorted(gff, prefix, fasta, out, out + ".idmap.tsv")
