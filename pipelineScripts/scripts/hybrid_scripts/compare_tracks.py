#!/usr/bin/env python3
"""Compare ARG recovery across the short / long / hybrid assembly tracks.

Consumes the output of `uniform_amr_compare.sbatch`, which runs one identical
ORF-calling + ARG-calling stack (prodigal -p meta -> AMRFinder+ -p / RGI -t protein
--local) over all three assemblies of the same sample. Because the calling stack is
held constant, differences reported here are attributable to the assemblies -- i.e.
to read type -- rather than to tool or database version.

Emits, under <outdir>:
  assembly_stats.tsv        per track x sample: contigs, bases, N50, ORFs
  arg_counts.tsv            per track x sample: AMRFinder+ / RGI call counts
  gene_matrix_<tool>.tsv    gene x (track,sample) presence matrix, per caller
  track_overlap.tsv         per sample: pairwise shared / unique gene counts + Jaccard
  track_unique_genes.tsv    genes recovered by exactly one track
  drug_class_by_track.tsv   RGI drug-class breakdown per track

Usage:
  python3 compare_tracks.py --amr-root /scratch/asorgen/hct_compare/amr \
                            --outdir  /scratch/asorgen/hct_compare/summary
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd

TRACKS = ["short", "long", "hybrid"]


# ---------------------------------------------------------------- assembly stats
def fasta_stats(path: Path) -> dict:
    """Contig count, total bases and N50 for a FASTA, in a single pass."""
    lengths, cur = [], 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)

    total = sum(lengths)
    n50, run = 0, 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run >= total / 2:
            n50 = L
            break
    return {"contigs": len(lengths), "bases": total, "n50": n50}


def count_orfs(faa: Path) -> int:
    if not faa.is_file():
        return 0
    with faa.open() as fh:
        return sum(1 for line in fh if line.startswith(">"))


# ------------------------------------------------------------------- ARG parsing
def read_amrfinder(path: Path) -> pd.DataFrame:
    """AMRFinder+ calls, restricted to actual AMR elements (not stress/virulence)."""
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["gene", "clazz"])
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        return pd.DataFrame(columns=["gene", "clazz"])
    if "Type" in df.columns:
        df = df[df["Type"].str.upper() == "AMR"]
    gene_col = "Element symbol" if "Element symbol" in df.columns else df.columns[1]
    class_col = "Class" if "Class" in df.columns else None
    out = pd.DataFrame({"gene": df[gene_col].str.strip()})
    out["clazz"] = df[class_col].str.strip() if class_col else ""
    return out[out["gene"] != ""]


def read_rgi(path: Path) -> pd.DataFrame:
    """RGI calls. Loose hits are dropped -- Strict/Perfect only."""
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["gene", "clazz"])
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        return pd.DataFrame(columns=["gene", "clazz"])
    if "Cut_Off" in df.columns:
        df = df[df["Cut_Off"].str.strip().isin(["Strict", "Perfect"])]
    gene_col = "Best_Hit_ARO" if "Best_Hit_ARO" in df.columns else df.columns[8]
    class_col = "Drug Class" if "Drug Class" in df.columns else None
    out = pd.DataFrame({"gene": df[gene_col].str.strip()})
    out["clazz"] = df[class_col].str.strip() if class_col else ""
    return out[out["gene"] != ""]


# -------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amr-root", required=True, type=Path,
                    help="root holding <track>/<sample>/ AMR outputs")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--samples", nargs="*", default=None,
                    help="restrict to these sample IDs (default: all found)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ---- discover completed (track, sample) pairs ----------------------------
    found: list[tuple[str, str]] = []
    for track in TRACKS:
        tdir = args.amr_root / track
        if not tdir.is_dir():
            continue
        for sdir in sorted(tdir.iterdir()):
            if not sdir.is_dir():
                continue
            if args.samples and sdir.name not in args.samples:
                continue
            if (sdir / "COMPLETE").exists():
                found.append((track, sdir.name))

    if not found:
        print(f"[ERROR] no completed track/sample pairs under {args.amr_root}",
              file=sys.stderr)
        return 1

    present_tracks = sorted({t for t, _ in found}, key=TRACKS.index)
    samples = sorted({s for _, s in found})
    print(f"[info] tracks: {', '.join(present_tracks)}")
    print(f"[info] samples: {', '.join(samples)}")
    incomplete = [(t, s) for t in present_tracks for s in samples
                  if (t, s) not in found]
    if incomplete:
        print(f"[warn] {len(incomplete)} track/sample pairs incomplete and excluded "
              f"from overlap: {incomplete}")

    # ---- per-(track, sample) stats and calls ---------------------------------
    stat_rows, count_rows, calls = [], [], {}
    for track, sample in found:
        d = args.amr_root / track / sample

        asm = d / f"{sample}_final_assembly.fasta"
        if not asm.is_file():  # short/long read their assembly from the pipeline tree
            asm = _pipeline_assembly(track, sample)
        st = fasta_stats(asm) if asm and asm.is_file() else {
            "contigs": 0, "bases": 0, "n50": 0}
        st.update(track=track, sample=sample, orfs=count_orfs(d / f"{sample}_genes.faa"))
        stat_rows.append(st)

        af = read_amrfinder(d / f"{sample}.amrfinder.tsv")
        rgi = read_rgi(d / f"{sample}.rgi.txt")
        calls[(track, sample, "amrfinder")] = af
        calls[(track, sample, "rgi")] = rgi
        count_rows.append({
            "track": track, "sample": sample,
            "amrfinder_calls": len(af), "amrfinder_genes": af["gene"].nunique(),
            "rgi_calls": len(rgi), "rgi_genes": rgi["gene"].nunique(),
        })

    stats = pd.DataFrame(stat_rows)[
        ["track", "sample", "contigs", "bases", "n50", "orfs"]
    ].sort_values(["sample", "track"])
    stats.to_csv(args.outdir / "assembly_stats.tsv", sep="\t", index=False)

    counts = pd.DataFrame(count_rows).sort_values(["sample", "track"])
    counts.to_csv(args.outdir / "arg_counts.tsv", sep="\t", index=False)

    # ---- gene x (track, sample) presence matrices ----------------------------
    for tool in ("amrfinder", "rgi"):
        recs = []
        for (track, sample, t), df in calls.items():
            if t != tool:
                continue
            for gene in df["gene"].unique():
                recs.append({"gene": gene, "col": f"{track}|{sample}", "present": 1})
        if recs:
            mat = (pd.DataFrame(recs)
                   .pivot_table(index="gene", columns="col", values="present",
                                fill_value=0, aggfunc="max")
                   .sort_index())
            mat.to_csv(args.outdir / f"gene_matrix_{tool}.tsv", sep="\t")

    # ---- pairwise track overlap, per sample ----------------------------------
    ov_rows, uniq_rows = [], []
    for tool in ("amrfinder", "rgi"):
        for sample in samples:
            sets = {
                tr: set(calls[(tr, sample, tool)]["gene"])
                for tr in present_tracks if (tr, sample) in found
            }
            for a, b in combinations(sorted(sets, key=TRACKS.index), 2):
                inter, union = sets[a] & sets[b], sets[a] | sets[b]
                ov_rows.append({
                    "tool": tool, "sample": sample, "track_a": a, "track_b": b,
                    "n_a": len(sets[a]), "n_b": len(sets[b]),
                    "shared": len(inter),
                    "only_a": len(sets[a] - sets[b]),
                    "only_b": len(sets[b] - sets[a]),
                    "jaccard": round(len(inter) / len(union), 4) if union else 0.0,
                })
            # genes seen in exactly one track for this sample
            if len(sets) > 1:
                for tr, genes in sets.items():
                    others = set().union(*(g for k, g in sets.items() if k != tr))
                    for gene in sorted(genes - others):
                        uniq_rows.append({"tool": tool, "sample": sample,
                                          "only_in_track": tr, "gene": gene})

    if ov_rows:
        pd.DataFrame(ov_rows).to_csv(args.outdir / "track_overlap.tsv",
                                     sep="\t", index=False)
    if uniq_rows:
        pd.DataFrame(uniq_rows).to_csv(args.outdir / "track_unique_genes.tsv",
                                       sep="\t", index=False)

    # ---- RGI drug-class breakdown per track ----------------------------------
    dc = []
    for (track, sample, tool), df in calls.items():
        if tool != "rgi" or df.empty:
            continue
        for classes in df["clazz"]:
            for c in (x.strip() for x in str(classes).split(";") if x.strip()):
                dc.append({"track": track, "sample": sample, "drug_class": c})
    if dc:
        (pd.DataFrame(dc)
         .groupby(["track", "drug_class"]).size()
         .reset_index(name="calls")
         .sort_values(["track", "calls"], ascending=[True, False])
         .to_csv(args.outdir / "drug_class_by_track.tsv", sep="\t", index=False))

    # ---- console summary ------------------------------------------------------
    print("\n=== assembly ===")
    print(stats.to_string(index=False))
    print("\n=== ARG calls ===")
    print(counts.to_string(index=False))
    if ov_rows:
        print("\n=== track overlap ===")
        print(pd.DataFrame(ov_rows).to_string(index=False))
    print(f"\n[done] wrote outputs to {args.outdir}")
    return 0


def _pipeline_assembly(track: str, sample: str) -> Path | None:
    """Locate the filtered assembly for tracks that read it from the pipeline tree."""
    study = Path("/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study")
    if track == "short":
        return (study / "HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_short"
                / "1.2_evaluation" / f"{sample}_final_assembly.fasta")
    if track == "long":
        return (study / "Nanopore_Resistome_Pipeline/results" / sample
                / "evaluate" / f"{sample}_final_assembly.fasta")
    return None


if __name__ == "__main__":
    sys.exit(main())
