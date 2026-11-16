#!/usr/bin/env python3
"""
ARG -> host MAG tables for one short-read cohort (bug-to-bug ARG transfer analysis).

For every sample:
  * every RGI hit on the assembly (5.5_AA_amr_assembly/<ID>/<ID>.rgi.txt, protein
    mode on Prodigal ORFs of the 1.2 final assembly; Strict/Perfect only, as RGI
    reports by default)
  * its contig -> refined MAG (3.2_refine_bins/<ID>/metawrap_<c>_<x>_bins/bin.N.fa,
    which keeps the original assembly contig names)
  * that MAG's reassembled version, quality (3.3 reassembled_bins.stats) and
    GTDB-Tk taxonomy (4.1 classify/gtdbtk.{bac120,ar53}.summary.tsv)

ARG names are RGI Best_Hit_ARO, i.e. the same CARD names RGI bwt uses (and that
module 4.14 uses), so genes join to the 4.14 results by name.

Unbinned ARG contigs are kept (MAG = NA): the fraction of ARG hits that land in a
MAG at all is part of the result.

Writes to <outdir>:
  <dataset>_mag_catalog.tsv   one row per MAG (incl. MAGs without ARGs)
  <dataset>_arg_hosts.tsv     one row per RGI hit on the assembly
  <dataset>_arg_host_build.log  per-sample counts and problems

Usage:
  build_arg_host_tables.py --dataset Heston_short --dataset-dir <datasetDir> \
      --outdir <PROCESSED_ROOT>/<dataset>_tables [--min-completion 70 --max-contam 10]
"""
import argparse
import csv
import glob
import os
import re
import sys

csv.field_size_limit(sys.maxsize)

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_gtdb(classification):
    """d__Bacteria;p__...;s__Genus species -> dict of rank -> name ('' if unassigned)."""
    out = {r: "" for r in RANKS}
    if not classification or classification.startswith("Unclassified"):
        return out
    for part, rank in zip(classification.split(";"), RANKS):
        out[rank] = part[3:] if len(part) > 3 else ""
    return out


def contig_of_orf(orf_id):
    """'NODE_3_length_299420_cov_48.952256_93 # 111434 # ...' -> 'NODE_3_length_299420_cov_48.952256'."""
    first = orf_id.split(" ")[0]
    return re.sub(r"_\d+$", "", first)


def contig_length(contig):
    m = re.search(r"_length_(\d+)_", contig)
    return int(m.group(1)) if m else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-completion", default="70")
    ap.add_argument("--max-contam", default="10")
    ap.add_argument("--limit", type=int, default=0, help="first N samples only (testing)")
    a = ap.parse_args()

    d = a.dataset_dir
    refined_tag = f"metawrap_{a.min_completion}_{a.max_contam}_bins"
    os.makedirs(a.outdir, exist_ok=True)

    ids = sorted(os.listdir(os.path.join(d, "5.5_AA_amr_assembly", "COMPLETE")))
    if a.limit:
        ids = ids[:a.limit]
    log = open(os.path.join(a.outdir, f"{a.dataset}_arg_host_build.log"), "w")
    log.write(f"{a.dataset}: {len(ids)} samples with 5.5 COMPLETE\n")
    log.write("SampleID\tn_MAGs\tn_ARG_hits\tn_ARG_hits_in_MAG\tnote\n")

    mag_cols = ["SampleID", "MAG", "MAG_version", "completeness", "contamination",
                "size_bp", "N50", "gtdb_classification"] + [f"gtdb_{r}" for r in RANKS] + \
               ["gtdb_fastani_reference", "gtdb_method"]
    arg_cols = ["SampleID", "ORF_ID", "contig", "contig_length", "Best_Hit_ARO", "ARO",
                "Cut_Off", "Best_Identities", "pct_length_of_ref", "Model_type",
                "AMR_Gene_Family", "Drug_Class", "Resistance_Mechanism",
                "MAG", "MAG_completeness", "MAG_contamination",
                "gtdb_family", "gtdb_genus", "gtdb_species"]

    mag_out = open(os.path.join(a.outdir, f"{a.dataset}_mag_catalog.tsv"), "w", newline="")
    arg_out = open(os.path.join(a.outdir, f"{a.dataset}_arg_hosts.tsv"), "w", newline="")
    mw = csv.DictWriter(mag_out, fieldnames=mag_cols, delimiter="\t", lineterminator="\n")
    aw = csv.DictWriter(arg_out, fieldnames=arg_cols, delimiter="\t", lineterminator="\n")
    mw.writeheader()
    aw.writeheader()

    tot_hits = tot_binned = tot_mags = 0
    for sid in ids:
        note = []

        # --- MAGs: reassembled version, quality, taxonomy ---
        mags = {}  # bin.N -> dict
        stats = os.path.join(d, "3.3_reassemble_bins", sid, "reassembled_bins.stats")
        if os.path.exists(stats):
            for r in read_tsv(stats):
                base = re.sub(r"\.(orig|strict|permissive)$", "", r["bin"])
                mags[base] = {"MAG_version": r["bin"], "completeness": r["completeness"],
                              "contamination": r["contamination"], "size_bp": r.get("size", ""),
                              "N50": r.get("N50", "")}
        tax = {}
        for dom in ("bac120", "ar53"):
            for f in (os.path.join(d, "4.1_classify_bins", sid, "classify", f"gtdbtk.{dom}.summary.tsv"),):
                if os.path.exists(f):
                    for r in read_tsv(f):
                        tax[r["user_genome"]] = r
        for base, m in mags.items():
            t = tax.get(m["MAG_version"])
            if t is None:
                note.append(f"no_gtdb:{m['MAG_version']}")
            cls = t.get("classification", "") if t else ""
            ranks = parse_gtdb(cls)
            m.update({"gtdb_classification": cls,
                      **{f"gtdb_{r}": ranks[r] for r in RANKS},
                      "gtdb_fastani_reference": (t or {}).get("closest_genome_reference",
                                                             (t or {}).get("fastani_reference", "")),
                      "gtdb_method": (t or {}).get("classification_method", "")})
            mw.writerow({"SampleID": sid, "MAG": base, **m})
        tot_mags += len(mags)

        # --- contig -> MAG, from the refined bins (original contig names) ---
        contig2mag = {}
        for fa in glob.glob(os.path.join(d, "3.2_refine_bins", sid, refined_tag, "bin.*.fa")):
            base = os.path.basename(fa)[:-3]
            with open(fa) as fh:
                for line in fh:
                    if line.startswith(">"):
                        contig2mag[line[1:].split()[0]] = base
        refined = set(contig2mag.values())
        if refined != set(mags):
            note.append(f"refined_vs_reassembled_mismatch:{len(refined)}vs{len(mags)}")

        # --- ARG hits on the assembly ---
        rgi = os.path.join(d, "5.5_AA_amr_assembly", sid, f"{sid}.rgi.txt")
        n = nb = 0
        if not os.path.exists(rgi):
            note.append("no_5.5_rgi")
        else:
            for r in read_tsv(rgi):
                contig = contig_of_orf(r["ORF_ID"])
                mag = contig2mag.get(contig, "")
                m = mags.get(mag, {})
                aw.writerow({
                    "SampleID": sid, "ORF_ID": r["ORF_ID"].split(" ")[0], "contig": contig,
                    "contig_length": contig_length(contig),
                    "Best_Hit_ARO": r["Best_Hit_ARO"], "ARO": r["ARO"], "Cut_Off": r["Cut_Off"],
                    "Best_Identities": r["Best_Identities"],
                    "pct_length_of_ref": r.get("Percentage Length of Reference Sequence", ""),
                    "Model_type": r["Model_type"], "AMR_Gene_Family": r["AMR Gene Family"],
                    "Drug_Class": r["Drug Class"], "Resistance_Mechanism": r["Resistance Mechanism"],
                    "MAG": mag or "NA",
                    "MAG_completeness": m.get("completeness", ""),
                    "MAG_contamination": m.get("contamination", ""),
                    "gtdb_family": m.get("gtdb_family", ""), "gtdb_genus": m.get("gtdb_genus", ""),
                    "gtdb_species": m.get("gtdb_species", "")})
                n += 1
                nb += bool(mag)
        tot_hits += n
        tot_binned += nb
        log.write(f"{sid}\t{len(mags)}\t{n}\t{nb}\t{';'.join(note)}\n")

    mag_out.close()
    arg_out.close()
    summary = (f"\nTOTAL: {len(ids)} samples, {tot_mags} MAGs, {tot_hits} ARG hits, "
               f"{tot_binned} in a MAG ({100 * tot_binned / max(tot_hits, 1):.1f}%)\n")
    log.write(summary)
    log.close()
    print(summary.strip())


if __name__ == "__main__":
    main()
