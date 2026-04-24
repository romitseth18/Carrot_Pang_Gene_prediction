#!/usr/bin/env python3

import sys
from collections import defaultdict
from intervaltree import IntervalTree

# ============================================================
# FASTA LOADER
# ============================================================
def load_fasta(fasta):
    seqs = {}
    name = None
    buf = []

    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

        if name:
            seqs[name] = "".join(buf).upper()

    return seqs


# ============================================================
# REVERSE COMPLEMENT
# ============================================================
def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


# ============================================================
# GFF PARSER
# ============================================================
def parse_gff(gff):
    features = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            p = line.strip().split("\t")
            if len(p) != 9:
                continue

            attrs = {}
            for a in p[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k] = v

            features.append({
                "seqid": p[0],
                "type": p[2],
                "start": int(p[3]),
                "end": int(p[4]),
                "strand": p[6],
                "attrs": attrs,
                "raw": p
            })

    return features


# ============================================================
# BUILD GENE INDEX (GFF2 backbone)
# ============================================================
def build_gene_index(gff2):
    idx = defaultdict(IntervalTree)

    for f in gff2:
        if f["type"] == "gene":
            gid = f["attrs"].get("ID")
            if gid:
                idx[f["seqid"]][f["start"]:f["end"]] = gid

    return idx


# ============================================================
# ASSIGN GENE (SAFE)
# ============================================================
def assign_gene(f, gene_index):
    if f["seqid"] not in gene_index:
        return None

    hits = gene_index[f["seqid"]].overlap(f["start"], f["end"])
    if not hits:
        return None

    best, score = None, 0
    for h in hits:
        s = min(f["end"], h.end) - max(f["start"], h.begin)
        if s > score:
            best, score = h.data, s

    return best


# ============================================================
# GROUP TRANSCRIPTS + FEATURES
# ============================================================
def group(gff1, gene_index):
    transcripts = defaultdict(list)
    children = defaultdict(list)

    for f in gff1:
        if f["type"] == "mRNA":
            gid = assign_gene(f, gene_index)
            tid = f["attrs"].get("ID")
            if gid and tid:
                transcripts[gid].append(f)

        elif f["type"] in ["CDS", "exon"]:
            parents = f["attrs"].get("Parent", "").split(",")
            for p in parents:
                if p:
                    children[p].append(f)

    return transcripts, children


# ============================================================
# BUILD CDS + PHASE (FRAME CORRECT)
# ============================================================
def build_cds_and_phase(features, genome):
    if not features:
        return None

    strand = features[0]["strand"]
    chrom = features[0]["seqid"]

    # Prefer CDS if available
    cds_feats = [f for f in features if f["type"] == "CDS"]
    if cds_feats:
        feats = cds_feats
    else:
        feats = features

    # Order
    if strand == "+":
        feats = sorted(feats, key=lambda x: x["start"])
    else:
        feats = sorted(feats, key=lambda x: x["start"], reverse=True)

    coords = []
    seq = ""

    for f in feats:
        frag = genome[chrom][f["start"]-1:f["end"]]
        seq += frag
        coords.append((f["start"], f["end"]))

    if strand == "-":
        seq = revcomp(seq)

    # Phase calculation
    phases = []
    frame = 0

    for start, end in coords:
        length = end - start + 1
        phases.append(frame)
        frame = (frame + length) % 3

    if strand == "-":
        coords.reverse()
        phases.reverse()

    return coords, phases


# ============================================================
# REBUILD GFF
# ============================================================
def rebuild(transcripts, children):
    out = []
    written_genes = set()

    for gene, txs in transcripts.items():

        for i, tx in enumerate(txs, 1):

            old_tx = tx["attrs"]["ID"]
            new_tx = f"{gene}_{i}"

            feats = children.get(old_tx, [])
            if not feats:
                continue

            cds_data = build_cds_and_phase(feats, genome)
            if not cds_data:
                continue

            coords, phases = cds_data
            seqid = feats[0]["seqid"]
            strand = feats[0]["strand"]

            # ---------------- GENE ----------------
            if gene not in written_genes:
                out.append("\t".join([
                    seqid, ".", "gene",
                    str(tx["start"]), str(tx["end"]),
                    ".", strand, ".", f"ID={gene}"
                ]))
                written_genes.add(gene)

            # ---------------- mRNA ----------------
            out.append("\t".join([
                seqid, ".", "mRNA",
                str(tx["start"]), str(tx["end"]),
                ".", strand, ".", f"ID={new_tx};Parent={gene}"
            ]))

            # ---------------- EXONS ----------------
            for e in feats:
                if e["type"] == "exon":
                    out.append("\t".join([
                        e["seqid"], ".", "exon",
                        str(e["start"]), str(e["end"]),
                        ".", e["strand"], ".",
                        f"Parent={new_tx}"
                    ]))

            # ---------------- CDS ----------------
            for (start, end), phase in zip(coords, phases):
                out.append("\t".join([
                    seqid, ".", "CDS",
                    str(start), str(end),
                    ".", strand, str(phase),
                    f"Parent={new_tx}"
                ]))

    return out


# ============================================================
# MAIN
# ============================================================
def main(gff1, gff2, fasta, out):

    print("[INFO] Loading genome...")
    global genome
    genome = load_fasta(fasta)

    print("[INFO] Parsing GFF1...")
    g1 = parse_gff(gff1)

    print("[INFO] Parsing GFF2...")
    g2 = parse_gff(gff2)

    print("[INFO] Building gene index...")
    gene_index = build_gene_index(g2)

    print("[INFO] Grouping transcripts...")
    transcripts, children = group(g1, gene_index)

    print("[DEBUG] genes in GFF2:", sum(1 for f in g2 if f["type"] == "gene"))
    print("[DEBUG] mapped transcripts:", len(transcripts))
    print("[DEBUG] transcript children:", len(children))

    print("[INFO] Rebuilding...")
    final = rebuild(transcripts, children)

    print("[INFO] Writing output...")
    with open(out, "w") as f:
        f.write("##gff-version 3\n")
        for l in final:
            f.write(l + "\n")

    print("[DONE]")


# ============================================================
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: script gff1 gff2 genome output")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
