#!/usr/bin/env python3
"""
GTDB genus -> NCBI genus lookup, from GTDB's own release metadata.

MAG hosts are named in GTDB (GTDB-Tk, release 226) and contig hosts in NCBI
(Kraken2), so the same organism often carries two genus names (e.g. GTDB
Thomasclavelia = NCBI Erysipelatoclostridium; GTDB Ruminococcus_B = NCBI
Mediterraneibacter). GTDB's metadata lists both taxonomies for every genome, so
for each GTDB genus we count the NCBI genera of its genomes.

Output columns (one row per GTDB genus, including suffixed ones like Blautia_A):
  gtdb_genus          GTDB genus as GTDB-Tk reports it
  gtdb_family         its GTDB family
  n_genomes           genomes in that GTDB genus with an NCBI genus
  ncbi_genus          majority NCBI genus
  majority_frac       share of those genomes with the majority NCBI genus
  ncbi_genera_ge10pct every NCBI genus covering >= 10% of the genomes ("|"-joined),
                      the set to use when asking "does a Kraken call match this MAG?"
  ncbi_genera_top     top 5 NCBI genera with counts, for inspection

With --kraken-taxonomy, NCBI genera are taken from the Kraken DB's own (older)
taxonomy snapshot by NCBI taxid, not from GTDB's current NCBI names, so the
lookup speaks the same names Kraken reports. "NCBI genus" below means that.

With --mag-catalogs and --kraken-tables it also reports how Kraken-vs-GTDB genus
agreement on binned ARG contigs changes when names go through the lookup.

Usage:
  build_gtdb_ncbi_genus_lookup.py --metadata bac120_metadata_r226.tsv.gz ar53_metadata_r226.tsv.gz \
      --out GTDB_r226_genus_to_NCBI_genus.tsv \
      [--mag-catalogs A_mag_catalog.tsv ... --kraken-tables A_arg_contig_kraken.tsv ...]
"""
import argparse
import collections
import csv
import gzip
import sys

csv.field_size_limit(sys.maxsize)


def rank(tax, prefix):
    for part in tax.split(";"):
        if part.startswith(prefix):
            return part[len(prefix):]
    return ""


class Taxonomy:
    """NCBI taxonomy snapshot (nodes/names/merged.dmp), e.g. the Kraken2 DB's own."""

    def __init__(self, taxdir):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        with open(f"{taxdir}/nodes.dmp") as fh:
            for line in fh:
                f = line.split("\t|\t")
                self.parent[f[0]] = f[1]
                self.rank[f[0]] = f[2]
        with open(f"{taxdir}/names.dmp") as fh:
            for line in fh:
                f = line.rstrip("\t|\n").split("\t|\t")
                if f[3] == "scientific name":
                    self.name[f[0]] = f[1]
        try:
            with open(f"{taxdir}/merged.dmp") as fh:
                for line in fh:
                    f = line.rstrip("\t|\n").split("\t|\t")
                    self.merged[f[0]] = f[1]
        except FileNotFoundError:
            pass

    def genus(self, taxid):
        t = self.merged.get(taxid, taxid)
        for _ in range(100):
            if t not in self.parent:
                return ""
            if self.rank[t] == "genus":
                return self.name.get(t, "")
            if self.parent[t] == t:
                return ""
            t = self.parent[t]
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--kraken-taxonomy",
                    help="taxonomy/ dir of the Kraken2 DB. Genomes are then mapped by NCBI taxid "
                         "through THAT snapshot, so genus names match what Kraken reports even "
                         "where NCBI has since renamed the genus (e.g. Erysipelatoclostridium -> "
                         "Thomasclavelia). Falls back to GTDB's ncbi_taxonomy string if the "
                         "taxid is newer than the snapshot.")
    ap.add_argument("--mag-catalogs", nargs="*", default=[])
    ap.add_argument("--kraken-tables", nargs="*", default=[])
    a = ap.parse_args()

    tax = Taxonomy(a.kraken_taxonomy) if a.kraken_taxonomy else None
    counts = collections.defaultdict(collections.Counter)
    family = {}
    via = collections.Counter()
    for path in a.metadata:
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gg = rank(r["gtdb_taxonomy"], "g__")
                if not gg:
                    continue
                family[gg] = rank(r["gtdb_taxonomy"], "f__")
                ng = ""
                if tax:
                    for col in ("ncbi_species_taxid", "ncbi_taxid"):
                        ng = tax.genus(r.get(col, "") or "")
                        if ng:
                            via["taxid"] += 1
                            break
                if not ng:
                    ng = rank(r["ncbi_taxonomy"], "g__")
                    via["ncbi_taxonomy_string" if ng else "none"] += 1
                if ng:
                    counts[gg][ng] += 1
    if tax:
        print(f"genome -> NCBI genus resolved via: {dict(via)}")
    lookup = {}
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["gtdb_genus", "gtdb_family", "n_genomes", "ncbi_genus", "majority_frac",
                    "ncbi_genera_ge10pct", "ncbi_genera_top"])
        for gg in sorted(family):
            c = counts.get(gg, collections.Counter())
            n = sum(c.values())
            if n:
                top, k = c.most_common(1)[0]
                ge10 = sorted(x for x, v in c.items() if v / n >= 0.10)
            else:
                top, k, ge10 = "", 0, []
            lookup[gg] = (top, set(ge10))
            w.writerow([gg, family[gg], n, top, f"{k / n:.3f}" if n else "", "|".join(ge10),
                        "; ".join(f"{x}({v})" for x, v in c.most_common(5))])
    print(f"{len(family)} GTDB genera written to {a.out} "
          f"({sum(1 for g in lookup if lookup[g][0])} with an NCBI genus)")

    # ---- optional: agreement on binned ARG contigs, before vs after the lookup ----
    for cat, kr in zip(a.mag_catalogs, a.kraken_tables):
        mags = {(r["SampleID"], r["MAG"]): r["gtdb_genus"]
                for r in csv.DictReader(open(cat), delimiter="\t")}
        raw = maj = anyset = tot = unmapped = 0
        left = collections.Counter()
        for r in csv.DictReader(open(kr), delimiter="\t"):
            if r["MAG"] == "NA" or not r["ncbi_genus"]:
                continue
            gg = mags.get((r["SampleID"], r["MAG"]), "")
            if not gg:
                continue
            tot += 1
            kg = r["ncbi_genus"]
            same_name = gg.split("_")[0] == kg
            raw += same_name
            top, ge10 = lookup.get(gg, ("", set()))
            if not top:
                unmapped += 1
            maj += same_name or kg == top
            ok = same_name or kg == top or kg in ge10
            anyset += ok
            if not ok:
                left[(gg, kg)] += 1
        print(f"\n{kr}\n  binned contigs with a GTDB and a Kraken genus: {tot}"
              f"\n  agree, suffix-stripped names only : {100 * raw / tot:.1f}%"
              f"\n  agree, via lookup (majority)      : {100 * maj / tot:.1f}%"
              f"\n  agree, via lookup (>=10% genera)  : {100 * anyset / tot:.1f}%"
              f"\n  MAG genus absent from lookup      : {unmapped}"
              f"\n  top remaining disagreements (GTDB -> Kraken): "
              + ", ".join(f"{g}->{k}({n})" for (g, k), n in left.most_common(8)))


if __name__ == "__main__":
    main()
