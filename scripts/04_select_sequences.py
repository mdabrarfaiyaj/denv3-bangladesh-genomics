"""
04_select_sequences.py
----------------------
Build the two analysis datasets from the Nextclade alignment.

Why two datasets: the whole-genome set has more signal per sequence but few taxa
(21 Bangladeshi genomes); the E-gene set has many more taxa (most older sequences
are E-gene fragments) but less signal. If both give the same answer, the answer is
much harder to argue with. They are not a main analysis plus a backup - they are two
independent tests.

Everything is driven by the Nextclade aligned FASTA, in which every sequence already
sits on reference coordinates (NC_001475.2). That means the E gene can be cut out by
position, and no extra alignment step is needed. MAFFT is kept as an optional
robustness check, not as a required step.

Inputs:
  results/metadata/metadata_clades.tsv
  results/nextclade/nextclade.aligned.fasta
  results/nextclade/dataset/genome_annotation.gff3   (for the E gene coordinates)

Outputs (in results/datasets/):
  wg.fasta          whole-genome alignment, selected taxa
  wg_dates.tsv      name, date      (TreeTime format, ambiguous dates allowed)
  wg_traits.tsv     name, country, region
  eg.fasta          E-gene alignment, selected taxa
  eg_dates.tsv
  eg_traits.tsv
  selection_log.tsv one row per sequence considered, with the reason it was kept or dropped

Run from the project folder:
  python scripts/04_select_sequences.py
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

METADATA = Path("results/metadata/metadata_clades.tsv")
ALIGNED = Path("results/nextclade/nextclade.aligned.fasta")
GFF = Path("results/nextclade/dataset/genome_annotation.gff3")
OUTDIR = Path("results/datasets")

# The focal lineage, its relatives used as context, and the ancestral lineage used
# to root the tree. 3III_A is included so the 2021 incursion can be analysed in the
# same framework, with its Indian relatives for comparison.
FOCAL = "3I_A.2"
CONTEXT = ["3I_A", "3I_A.1"]
OUTGROUP = ["3I"]
INCURSION = "3III_A"

# ---------------------------------------------------------------------------
# Decisions from the curation audit (08) and the recombination work (09, 10).
# Each is recorded here with its reason so the dataset is self-documenting.
# See claude/project1_denv3_decisions_log.md.
# ---------------------------------------------------------------------------

# D1. Confirmed phylogenetic incongruence: PQ815901.1 groups with genotype III
# in three blocks (1501-2250, 3751-4250, 4501-5250) and with genotype I in the
# rest, both at SH-aLRT 100 / UFBoot 100, while no control sequence moves. One
# tree cannot represent a genome with two histories, so it is excluded.
EXCLUDE = {
    "PQ815901.1": "probable recombinant / chimeric assembly (confirmed incongruence)",
}

# D6. Rescued from the frameshift filter. Its two frameshifts are both in E
# (E:192-193, E:394) and every other gene is at coverage 1.0, so E is masked for
# this sequence only. As the earliest 3I_A.2 whole genome (China, 2015-12-22) it
# is the best available handle on the origin date, and dropping it on an
# unexamined QC flag was a curation error.
RESCUE = {
    "MN018377.1": ["E"],
}

# D5. Month-end dates from a single submitter are end-of-collection-period
# placeholders, not collection days. Demoted to month precision rather than
# being trusted as exact or thrown away.
DATE_OVERRIDE = {
    **{f"PZ0174{n}.1": ("2025-09", "month") for n in range(40, 45)},
    **{f"PZ0174{n}.1": ("2024-10", "month") for n in range(35, 40)},
}

# Coarse regions for the mugration step. Country-level mugration with many
# single-sequence countries is noisy, so both levels are written out and the
# analysis reports both.
REGION = {
    "Bangladesh": "Bangladesh",
    "India": "South Asia (other)",
    "Sri Lanka": "South Asia (other)",
    "Pakistan": "South Asia (other)",
    "Nepal": "South Asia (other)",
    "Malaysia": "Southeast Asia",
    "Indonesia": "Southeast Asia",
    "Singapore": "Southeast Asia",
    "Thailand": "Southeast Asia",
    "Philippines": "Southeast Asia",
    "Vietnam": "Southeast Asia",
    "Cambodia": "Southeast Asia",
    "Myanmar": "Southeast Asia",
    "Laos": "Southeast Asia",
    "East Timor": "Southeast Asia",
    "Brunei": "Southeast Asia",
    "China": "East Asia",
    "Taiwan": "East Asia",
    "Japan": "East Asia",
    "South Korea": "East Asia",
    "Hong Kong": "East Asia",
}


def region_of(country):
    if pd.isna(country):
        return "unknown"
    return REGION.get(country, "elsewhere")


def treetime_date(date, precision):
    """TreeTime accepts ambiguous dates as YYYY-XX-XX. Never invent precision."""
    if pd.isna(date) or precision == "missing":
        return None
    if precision == "day":
        return date
    if precision == "month":
        return f"{date}-XX"
    return f"{date}-XX-XX"


def find_e_gene(gff_path):
    """Read the E gene coordinates from the dataset annotation, rather than hardcoding."""
    if not gff_path.exists():
        return None, []
    feats = []
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            start, end, attrs = int(f[3]), int(f[4]), f[8]
            m = re.search(r"(?:gene_name|Name|gene|ID)=([^;]+)", attrs)
            name = m.group(1) if m else "?"
            feats.append((name, f[2], start, end))
    e = [x for x in feats if re.fullmatch(r"E|envelope|E_gene", x[0], flags=re.I)]
    return (e[0] if e else None), feats


def coverage(seq):
    """Fraction of the slice that is a real base."""
    if not seq:
        return 0.0
    good = sum(1 for c in seq.upper() if c in "ACGT")
    return good / len(seq)


def read_gene_coords(gff_path):
    """{gene: (nt_start, nt_end)} from the dataset annotation, for masking."""
    genes = {}
    if not gff_path.exists():
        return genes
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = re.search(r"(?:gene_name|Name|gene|ID)=([^;]+)", f[8])
            if m:
                genes.setdefault(m.group(1), (int(f[3]), int(f[4])))
    return genes


def mask_genes(seq, genes, gene_coords):
    """Replace the named genes with N for one sequence, keeping alignment length."""
    s = list(seq)
    for g in genes:
        gc = gene_coords.get(g)
        if not gc:
            continue
        for i in range(gc[0] - 1, min(gc[1], len(s))):
            s[i] = "N"
    return "".join(s)


def subsample(df, per_group, seed):
    """Cap the number of context sequences per country and year, reproducibly."""
    if df.empty:
        return df
    return (df.sort_values("accession")
              .groupby(["country", "year"], dropna=False, group_keys=False)
              .apply(lambda g: g.sample(min(len(g), per_group), random_state=seed)))


def main():
    # declared up front: the argparse defaults below READ these module-level
    # values, and Python forbids a global declaration after first use
    global OUTDIR, OUTGROUP

    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="subsampling seed; rerun with other seeds as a robustness check")
    ap.add_argument("--context-per-group", type=int, default=2,
                    help="max context sequences per country per year")
    ap.add_argument("--min-cov-wg", type=float, default=0.80,
                    help="minimum fraction of real bases across the genome")
    ap.add_argument("--min-cov-eg", type=float, default=0.70,
                    help="minimum fraction of real bases across the E gene")
    ap.add_argument("--outdir", default=str(OUTDIR),
                    help="where to write the datasets; use a separate directory "
                         "for each robustness run so nothing is overwritten")
    ap.add_argument("--outgroup-lineages", default=",".join(OUTGROUP),
                    help="comma-separated lineages to use as the outgroup. The "
                         "default '3I' is a PARAPHYLETIC parent label and cannot "
                         "root the tree; pass a sister genotype such as '3II_B' "
                         "for a genuine outgroup rooting.")
    ap.add_argument("--outgroup-max", type=int, default=0,
                    help="cap on outgroup sequences (0 = no cap). Keep this small "
                         "when using a sister genotype: the outgroup is there to "
                         "root the tree, not to change the analysis.")
    args = ap.parse_args()

    OUTDIR = Path(args.outdir)
    OUTGROUP = [x.strip() for x in args.outgroup_lineages.split(",") if x.strip()]
    print(f"Output directory : {OUTDIR}")
    print(f"Outgroup lineages: {OUTGROUP}"
          f"{f' (capped at {args.outgroup_max})' if args.outgroup_max else ''}")
    if OUTGROUP == ["3I"]:
        print("  NOTE: '3I' is a paraphyletic parent label, so IQ-TREE cannot "
              "root on it.\n        The tree will be clock-rooted instead. Pass "
              "--outgroup-lineages 3II_B for a true outgroup.")

    for p in (METADATA, ALIGNED):
        if not p.exists():
            sys.exit(f"Input not found: {p}")

    md = pd.read_csv(METADATA, sep="\t", low_memory=False)

    e_feat, all_feats = find_e_gene(GFF)
    print("== Annotated features in the Nextclade dataset ==")
    for name, kind, s, en in all_feats:
        print(f"  {name:12s} {kind:12s} {s:6d}-{en:6d}")
    if e_feat is None:
        sys.exit("\nCould not find an E gene feature in the annotation. "
                 "Inspect the list above and set the coordinates explicitly.")
    _, _, e_start, e_end = e_feat
    print(f"\nE gene taken from the annotation: {e_start}-{e_end} "
          f"({e_end - e_start + 1} bases)")

    # ---- who is eligible, on metadata grounds alone ----
    # D5: demote month-end placeholder dates to month precision before any
    # TreeTime-format date is built from them.
    applied = []
    for acc, (newdate, newprec) in DATE_OVERRIDE.items():
        hit = md["accession"] == acc
        if hit.any():
            old = md.loc[hit, "date"].iloc[0]
            md.loc[hit, "date"] = newdate
            md.loc[hit, "date_precision"] = newprec
            applied.append(f"{acc}: {old} -> {newdate} ({newprec})")
    if applied:
        print(f"\nDate precision demoted for {len(applied)} sequences "
              f"(month-end placeholders, see decisions log D5):")
        for line in applied:
            print(f"  {line}")

    keep_clades = [FOCAL, INCURSION] + CONTEXT + OUTGROUP
    md["ttdate"] = [treetime_date(d, p) for d, p in zip(md["date"], md["date_precision"])]
    md["region"] = md["country"].map(region_of)

    elig = md[md["clade"].isin(keep_clades) & md["ttdate"].notna() & md["country"].notna()].copy()
    print(f"\nSequences in the focal + context + outgroup lineages, dated and with a country: {len(elig)}")

    # D1: exclude sequences disqualified on the evidence, with the reason printed
    # so the dataset explains itself.
    for acc, why in EXCLUDE.items():
        hit = elig["accession"] == acc
        if hit.any():
            print(f"\nEXCLUDED {acc}: {why}")
            elig = elig[~hit]

    # Frameshifts or premature stop codons in a viral coding sequence almost always
    # mean an assembly error, so those sequences are dropped rather than masked -
    # except where D6 rescues one and masks only the affected gene.
    for col in ("qc.frameShifts.status", "qc.stopCodons.status"):
        if col in elig.columns:
            bad = elig[col].eq("bad") & ~elig["accession"].isin(RESCUE)
            rescued = elig[col].eq("bad") & elig["accession"].isin(RESCUE)
            if bad.any():
                print(f"Dropping {bad.sum()} sequences flagged bad for {col}")
            if rescued.any():
                for acc in elig.loc[rescued, "accession"]:
                    print(f"RESCUED {acc} despite {col}=bad; masking "
                          f"{', '.join(RESCUE[acc])} (see decisions log D6)")
            elig = elig[~bad]

    want = set(elig["accession"])

    # ---- read the alignment once, keeping only what might be used ----
    print("\nReading the alignment ...")
    gene_coords = read_gene_coords(GFF)
    wg_seq, eg_seq = {}, {}
    for header, seq in SimpleFastaParser(ALIGNED.open()):
        acc = header.split()[0]
        if acc not in want:
            continue
        if acc in RESCUE:
            seq = mask_genes(seq, RESCUE[acc], gene_coords)
        wg_seq[acc] = seq
        eg_seq[acc] = seq[e_start - 1:e_end]
    print(f"Found {len(wg_seq)} of {len(want)} eligible sequences in the alignment")

    elig["cov_wg"] = elig["accession"].map(lambda a: coverage(wg_seq.get(a, "")))
    elig["cov_eg"] = elig["accession"].map(lambda a: coverage(eg_seq.get(a, "")))

    print("\nCoverage of the Bangladeshi sequences (how much real sequence there is):")
    bd = elig[elig["country"] == "Bangladesh"]
    print(bd.groupby("year")[["cov_wg", "cov_eg"]].mean().round(3).to_string())

    # ---- dataset A: whole genome ----
    a_pass = elig[elig["cov_wg"] >= args.min_cov_wg].copy()
    a_focal = a_pass[a_pass["clade"].isin([FOCAL, INCURSION])]
    a_ctx = subsample(a_pass[a_pass["clade"].isin(CONTEXT)], args.context_per_group, args.seed)
    a_out = a_pass[a_pass["clade"].isin(OUTGROUP)]
    if args.outgroup_max and len(a_out) > args.outgroup_max:
        a_out = a_out.sort_values("accession").sample(args.outgroup_max,
                                                      random_state=args.seed)
    A = pd.concat([a_focal, a_ctx, a_out]).drop_duplicates("accession")

    # ---- dataset B: E gene ----
    b_pass = elig[elig["cov_eg"] >= args.min_cov_eg].copy()
    b_focal = b_pass[b_pass["clade"].isin([FOCAL, INCURSION])]
    b_ctx = subsample(b_pass[b_pass["clade"].isin(CONTEXT)], args.context_per_group, args.seed)
    b_out = b_pass[b_pass["clade"].isin(OUTGROUP)]
    if args.outgroup_max and len(b_out) > args.outgroup_max:
        b_out = b_out.sort_values("accession").sample(args.outgroup_max,
                                                      random_state=args.seed)
    B = pd.concat([b_focal, b_ctx, b_out]).drop_duplicates("accession")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    def write_set(df, seqs, stem):
        with (OUTDIR / f"{stem}.fasta").open("w") as fh:
            for acc in df["accession"]:
                fh.write(f">{acc}\n{seqs[acc]}\n")
        df[["accession", "ttdate"]].rename(
            columns={"accession": "name", "ttdate": "date"}
        ).to_csv(OUTDIR / f"{stem}_dates.tsv", sep="\t", index=False)
        t = df[["accession", "country", "region", "clade"]].rename(
            columns={"accession": "name"}).copy()
        # mark the outgroup here, so downstream scripts never have to know which
        # lineage was chosen for it
        t["is_outgroup"] = t["clade"].isin(OUTGROUP).map({True: "yes", False: "no"})
        t.to_csv(OUTDIR / f"{stem}_traits.tsv", sep="\t", index=False)

    write_set(A, wg_seq, "wg")
    write_set(B, eg_seq, "eg")

    # ---- log every decision ----
    log = elig[["accession", "country", "region", "year", "date_precision", "clade",
                "length", "cov_wg", "cov_eg", "qc.overallStatus"]].copy()
    log["in_wg"] = log["accession"].isin(set(A["accession"]))
    log["in_eg"] = log["accession"].isin(set(B["accession"]))
    log.to_csv(OUTDIR / "selection_log.tsv", sep="\t", index=False)

    # ---- report ----
    for name, D in (("WHOLE GENOME", A), ("E GENE", B)):
        print(f"\n================ {name}: {len(D)} sequences ================")
        print("\nBy lineage:")
        print(D["clade"].value_counts().to_string())
        print("\nBangladesh by year:")
        bdd = D[D["country"] == "Bangladesh"]
        if len(bdd):
            print(pd.crosstab(bdd["year"].astype(int), bdd["clade"]).to_string())
        print(f"\nBangladesh: {len(bdd)}   elsewhere: {len(D) - len(bdd)}")
        print("Top context countries:")
        print(D[D["country"] != "Bangladesh"]["country"].value_counts().head(10).to_string())
        print("Date precision:")
        print(D["date_precision"].value_counts().to_string())

    print(f"\nWritten to {OUTDIR}/  (seed={args.seed}, context cap={args.context_per_group})")
    print("Check that every 2025 Bangladeshi sequence you expect is present before moving on.")


if __name__ == "__main__":
    main()
