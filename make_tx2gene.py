from pathlib import Path
import re

gtf_path = Path("Saccharomyces_cerevisiae.R64-1-1.115.gtf")
out_path = Path("tx2gene.tsv")

attr_re = re.compile(r'(\S+) "([^"]+)"')

with gtf_path.open() as fin, out_path.open("w") as fout:
    fout.write("transcript_id\tgene_id\n")
    seen = set()

    for line in fin:
        if line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        feature = cols[2]
        attrs = cols[8]

        # transcript 行を優先、なければ exon でも抽出可
        if feature not in {"transcript", "mRNA", "exon", "CDS"}:
            continue

        attr_map = dict(attr_re.findall(attrs))
        tx = attr_map.get("transcript_id")
        gene = attr_map.get("gene_id")

        if tx and gene and (tx, gene) not in seen:
            fout.write(f"{tx}\t{gene}\n")
            seen.add((tx, gene))

print(f"written: {out_path}")
