"""
05_mutation_profile.py
----------------------
What is this virus, at the protein level?

Nextclade already called every amino acid substitution against the DENV-3 reference
(NC_001475.2). This script mines those calls to answer three questions:

  1. Which substitutions define the Bangladeshi 3I_A.2 lineage?
  2. Which substitutions separate the 2025 sequences from the 2017-2022 ones?
     (i.e. did the persisting lineage change, and where?)
  3. Which are Bangladesh-specific, and which are shared with the Southeast Asian
     relatives the lineage descends from?

The coverage problem, handled exactly (this was a bug in the first version):
A substitution cannot be called absent at a position that was never sequenced.
Nextclade's `unknownAaRanges` is NOT sufficient for this: it reports unknown
residues INSIDE the aligned region only, and says nothing about genes lying
outside alignmentStart-alignmentEnd. PZ017440.1 shows the problem - its
cdsCoverage says NS1:0.088, but unknownAaRanges lists only C:87, so NS1 looks
fully sequenced when 91% of it is absent.

So a sequence counts toward the denominator at gene G residue p only if
  (a) the codon for p lies inside [alignmentStart, alignmentEnd], computed from
      the gene's nucleotide coordinates in the dataset annotation, AND
  (b) p is not inside an unknownAaRanges interval for G.
Both the count and the denominator are reported for every group. A frequency with
a small denominator means "we hardly looked", not "it is rare".

Input : results/nextclade/nextclade.tsv
        results/metadata/metadata_clades.tsv
        results/nextclade/dataset/genome_annotation.gff3   (gene coordinates)
Output: results/mutations/substitution_table.tsv   every substitution, by group
        results/mutations/bd_by_year.tsv           Bangladesh frequency per year
        results/mutations/swept.tsv                substitutions that rose to fixation
        results/mutations/lineage_defining.tsv     the Bangladesh 3I_A.2 signature
        results/mutations/e_protein_domains.tsv    E substitutions with domain annotation

Run from the project folder:
  python scripts/05_mutation_profile.py
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

NEXTCLADE = Path("results/nextclade/nextclade.tsv")
METADATA = Path("results/metadata/metadata_clades.tsv")
GFF = Path("results/nextclade/dataset/genome_annotation.gff3")
OUTDIR = Path("results/mutations")

# Groups compared. Order matters only for the printout.
GROUPS = {
    "BD_2025":      lambda d: (d.country == "Bangladesh") & (d.year == 2025),
    "BD_2021_2024": lambda d: (d.country == "Bangladesh") & d.year.between(2021, 2024),
    "BD_2017_2019": lambda d: (d.country == "Bangladesh") & d.year.between(2017, 2019),
    "3I_A.2_other": lambda d: (d.country != "Bangladesh") & (d.clade == "3I_A.2"),
    "3I_A_context": lambda d: d.clade.isin(["3I_A", "3I_A.1"]),
}

# Envelope protein domain boundaries, standard flavivirus E numbering (residues).
# DI is discontinuous; DII carries the fusion loop; DIII is the receptor-binding
# surface and the target of the most potently neutralising antibodies.
E_DOMAINS = [
    (1, 52, "domain I"),
    (53, 131, "domain II"),
    (132, 193, "domain I"),
    (194, 279, "domain II"),
    (280, 296, "domain I"),
    (297, 394, "domain III"),
    (395, 495, "stem / anchor"),
]
E_FEATURES = [
    (98, 110, "fusion loop"),
    (67, 67, "N-glycosylation site (N67)"),
    (153, 153, "N-glycosylation site (N153)"),
    (62, 80, "glycan loop"),
]


def e_annotation(pos):
    dom = next((n for a, b, n in E_DOMAINS if a <= pos <= b), "outside annotated range")
    feats = [n for a, b, n in E_FEATURES if a <= pos <= b]
    return dom, "; ".join(feats)


def parse_subs(cell):
    """'E:T329I,NS5:K12R' -> [('E', 329, 'T329I'), ('NS5', 12, 'K12R')]"""
    if not isinstance(cell, str) or not cell.strip():
        return []
    out = []
    for tok in cell.split(","):
        tok = tok.strip()
        if ":" not in tok:
            continue
        gene, change = tok.split(":", 1)
        m = re.match(r"^([A-Z*\-])(\d+)([A-Z*\-])$", change)
        if m:
            out.append((gene, int(m.group(2)), change))
    return out


def parse_unknown_ranges(cell):
    """
    Nextclade's unknownAaRanges looks like 'E:21-211,C:87' - residues not read.
    Returns {gene: [(start, end), ...]} so we know where we were blind.
    Entries can be a range ('E:21-211') or a single residue ('C:87').
    """
    cov = defaultdict(list)
    if not isinstance(cell, str) or not cell.strip():
        return cov
    for tok in cell.split(","):
        tok = tok.strip()
        m = re.match(r"^([A-Za-z0-9_]+):(\d+)-(\d+)$", tok)
        if m:
            cov[m.group(1)].append((int(m.group(2)), int(m.group(3))))
            continue
        m = re.match(r"^([A-Za-z0-9_]+):(\d+)$", tok)
        if m:
            p = int(m.group(2))
            cov[m.group(1)].append((p, p))
    return cov


def read_gene_coords(gff_path):
    """{gene: (nt_start, nt_end)} from the dataset annotation."""
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


def covered(gene, pos, unknown, aln_start, aln_end, gene_coords):
    """
    True only if this sequence really has data at gene residue `pos`.

    Two independent ways it can be missing:
      1. the residue's codon lies outside the aligned region entirely
         (unknownAaRanges does NOT report this - the bug this fixes)
      2. the residue is an N inside the aligned region (unknownAaRanges does
         report this)
    """
    gc = gene_coords.get(gene)
    if gc is None:
        return False                      # unknown gene: refuse to count it
    if pd.isna(aln_start) or pd.isna(aln_end):
        return False
    codon_start = gc[0] + (pos - 1) * 3
    codon_end = codon_start + 2
    if codon_start < aln_start or codon_end > aln_end:
        return False
    for a, b in unknown.get(gene, []):
        if a <= pos <= b:
            return False
    return True


def main():
    for p in (NEXTCLADE, METADATA):
        if not p.exists():
            sys.exit(f"Input not found: {p}")

    nc = pd.read_csv(NEXTCLADE, sep="\t", low_memory=False)
    md = pd.read_csv(METADATA, sep="\t", low_memory=False)

    sub_col = next((c for c in ["aaSubstitutions", "aaSubstitutions.1"] if c in nc.columns), None)
    if sub_col is None:
        print("Columns available in nextclade.tsv:")
        print("\n".join(f"  {c}" for c in nc.columns))
        sys.exit("\nNo aaSubstitutions column found. Paste the list above.")
    unk_col = next((c for c in ["unknownAaRanges", "aaUnknowns", "unknownAaRanges.1"]
                    if c in nc.columns), None)
    for need in ("alignmentStart", "alignmentEnd"):
        if need not in nc.columns:
            sys.exit(f"nextclade.tsv has no {need} column; cannot compute honest "
                     f"coverage denominators.")

    gene_coords = read_gene_coords(GFF)
    if not gene_coords:
        sys.exit(f"Could not read gene coordinates from {GFF}.")
    print("Gene coordinates used for coverage (from the dataset annotation):")
    for g, (a, b) in sorted(gene_coords.items(), key=lambda kv: kv[1]):
        print(f"  {g:6s} {a:6d}-{b:6d}  ({(b - a + 1) // 3} residues)")

    nc["accession"] = nc["seqName"].astype(str).str.split().str[0]

    # metadata_clades.tsv already carries alignmentStart/alignmentEnd from step 2b,
    # so an ordinary merge silently produces alignmentStart_x / _y. Rename the
    # Nextclade columns explicitly instead of trusting merge suffixes.
    rename = {sub_col: "aaSubs_nc",
              "alignmentStart": "alnStart_nc",
              "alignmentEnd": "alnEnd_nc"}
    if unk_col:
        rename[unk_col] = "unkAa_nc"
    ncsub = nc[["accession"] + list(rename)].rename(columns=rename)
    d = md.drop(columns=[c for c in rename.values() if c in md.columns],
                errors="ignore").merge(ncsub, on="accession", how="inner")
    d["year"] = pd.to_numeric(d["year"], errors="coerce")

    # Only genotype I matters here; 3III_A is a separate story handled elsewhere.
    d = d[d["clade"].isin(["3I_A.2", "3I_A", "3I_A.1", "3I"])].copy()
    print(f"Sequences considered: {len(d)}")

    # ---- tally substitutions per group, with honest denominators ----
    counts = defaultdict(lambda: defaultdict(int))   # sub -> group -> n with it
    denom = defaultdict(lambda: defaultdict(int))    # sub -> group -> n that covered it
    gene_of, pos_of = {}, {}

    group_of = {}
    for gname, pred in GROUPS.items():
        for acc in d[pred(d)]["accession"]:
            group_of[acc] = gname

    # first pass: learn which substitutions exist at all, and where
    all_subs = set()
    parsed, unknowns, aln = {}, {}, {}
    for row in d.itertuples(index=False):
        acc = row.accession
        subs = parse_subs(row.aaSubs_nc)
        parsed[acc] = subs
        unknowns[acc] = (parse_unknown_ranges(row.unkAa_nc) if unk_col
                         else defaultdict(list))
        aln[acc] = (row.alnStart_nc, row.alnEnd_nc)
        for gene, pos, change in subs:
            key = f"{gene}:{change}"
            all_subs.add(key)
            gene_of[key], pos_of[key] = gene, pos

    def is_covered(acc, key):
        a, b = aln.get(acc, (None, None))
        return covered(gene_of[key], pos_of[key], unknowns.get(acc, {}), a, b, gene_coords)

    # second pass: count, skipping sequences blind at that position
    by_year = defaultdict(lambda: defaultdict(lambda: [0, 0]))   # sub -> year -> [n, of]
    bd_accs = set(d[d["country"] == "Bangladesh"]["accession"])
    year_of = dict(zip(d["accession"], d["year"]))

    for acc, g in group_of.items():
        subs_here = {f"{gene}:{ch}" for gene, _, ch in parsed.get(acc, [])}
        for key in all_subs:
            if not is_covered(acc, key):
                continue
            denom[key][g] += 1
            if key in subs_here:
                counts[key][g] += 1

    # Bangladesh frequency per year: the three-period grouping hides sweeps
    for acc in bd_accs:
        yr = year_of.get(acc)
        if pd.isna(yr):
            continue
        yr = int(yr)
        subs_here = {f"{gene}:{ch}" for gene, _, ch in parsed.get(acc, [])}
        for key in all_subs:
            if not is_covered(acc, key):
                continue
            by_year[key][yr][1] += 1
            if key in subs_here:
                by_year[key][yr][0] += 1

    rows = []
    for key in sorted(all_subs, key=lambda k: (gene_of[k], pos_of[k])):
        r = {"substitution": key, "gene": gene_of[key], "position": pos_of[key]}
        for g in GROUPS:
            n, m = counts[key].get(g, 0), denom[key].get(g, 0)
            r[f"{g}_n"] = n
            r[f"{g}_of"] = m
            r[f"{g}_freq"] = round(n / m, 3) if m else None
        rows.append(r)
    tab = pd.DataFrame(rows)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUTDIR / "substitution_table.tsv", sep="\t", index=False)
    print(f"All substitutions -> {OUTDIR / 'substitution_table.tsv'}  ({len(tab)} rows)")

    for g in GROUPS:
        n = sum(1 for a, gg in group_of.items() if gg == g)
        print(f"  group {g:14s} {n:4d} sequences")

    # ---- 1a. Bangladesh frequency per year ----
    years = sorted({y for v in by_year.values() for y in v})
    yr_rows = []
    for key in sorted(all_subs, key=lambda k: (gene_of[k], pos_of[k])):
        r = {"substitution": key, "gene": gene_of[key], "position": pos_of[key]}
        for y in years:
            n, m = by_year[key].get(y, [0, 0])
            r[f"{y}_n"] = n
            r[f"{y}_of"] = m
            r[f"{y}"] = round(n / m, 3) if m else None
        yr_rows.append(r)
    yrtab = pd.DataFrame(yr_rows)
    yrtab.to_csv(OUTDIR / "bd_by_year.tsv", sep="\t", index=False)
    print(f"\nBangladesh frequency per year -> {OUTDIR / 'bd_by_year.tsv'}")
    print(f"Years present: {', '.join(map(str, years))}")

    # ---- 1b. substitutions that swept to fixation inside Bangladesh ----
    # Low early, high late, with enough coverage at both ends. This is what a
    # real selective or founder sweep looks like, and the three-period grouping
    # in the first version of this script hid it.
    early = [y for y in years if y <= 2018]
    late = [y for y in years if y >= 2021]
    if early and late:
        e_n = yrtab[[f"{y}_n" for y in early]].sum(axis=1)
        e_m = yrtab[[f"{y}_of" for y in early]].sum(axis=1)
        l_n = yrtab[[f"{y}_n" for y in late]].sum(axis=1)
        l_m = yrtab[[f"{y}_of" for y in late]].sum(axis=1)
        sw = yrtab.copy()
        sw["early_freq"] = (e_n / e_m).round(3)
        sw["early_of"] = e_m
        sw["late_freq"] = (l_n / l_m).round(3)
        sw["late_of"] = l_m
        swept = sw[(sw["early_of"] >= 10) & (sw["late_of"] >= 10) &
                   (sw["early_freq"] <= 0.25) & (sw["late_freq"] >= 0.90)]
        swept.to_csv(OUTDIR / "swept.tsv", sep="\t", index=False)
        print(f"\n== Rose from <=0.25 (by 2018) to >=0.90 (2021 on) in Bangladesh: "
              f"{len(swept)} ==")
        if len(swept):
            print(swept[["substitution", "gene", "position",
                         "early_freq", "early_of", "late_freq", "late_of"]]
                  .to_string(index=False))

    # ---- 1c. the Bangladesh lineage signature ----
    # Fixed in Bangladeshi sequences and clearly rarer in the Southeast Asian
    # relatives. Two thresholds are reported because the strict one returned
    # nothing: within-country heterogeneity in 2017-2019 is real, not noise.
    bd_cols = ["BD_2017_2019", "BD_2021_2024", "BD_2025"]
    sig = tab.copy()
    sig["bd_min_freq"] = sig[[f"{c}_freq" for c in bd_cols]].min(axis=1, skipna=True)
    sig["bd_recent_freq"] = sig[[f"{c}_freq" for c in ["BD_2021_2024", "BD_2025"]]].min(axis=1, skipna=True)
    sig["bd_covered"] = sig[[f"{c}_of" for c in bd_cols]].sum(axis=1)

    strict = sig[(sig["bd_min_freq"] >= 0.90) & (sig["bd_covered"] >= 10) &
                 (sig["3I_A_context_freq"].fillna(0) <= 0.10)]
    recent = sig[(sig["bd_recent_freq"] >= 0.90) & (sig["bd_covered"] >= 10) &
                 (sig["3I_A_context_freq"].fillna(0) <= 0.10)]
    recent.to_csv(OUTDIR / "lineage_defining.tsv", sep="\t", index=False)

    cols = ["substitution", "gene", "position", "bd_min_freq", "bd_recent_freq",
            "3I_A.2_other_freq", "3I_A_context_freq"]
    print(f"\n== Signature, strict (fixed in every Bangladeshi period): {len(strict)} ==")
    if len(strict):
        print(strict[cols].to_string(index=False))
    print(f"\n== Signature, 2021 onwards (fixed in the current lineage): {len(recent)} ==")
    if len(recent):
        print(recent[cols].to_string(index=False))
    print("\nA substitution high in '3I_A.2_other' was INHERITED when the lineage")
    print("arrived (~2015-16); one near zero there arose in Bangladesh.")

    # ---- 2. what changed by 2025 ----
    # Present in 2025, rare or absent in 2017-2019. Coverage denominators matter
    # most here, because the 2025 genomes are the incomplete ones.
    chg = tab[(tab["BD_2025_freq"].fillna(0) >= 0.60) &
              (tab["BD_2025_of"] >= 3) &
              (tab["BD_2017_2019_freq"].fillna(0) <= 0.20) &
              (tab["BD_2017_2019_of"] >= 5)]
    print(f"\n== Substitutions gained by 2025 (vs 2017-2019): {len(chg)} ==")
    print("Read the _of columns: a low denominator means thin evidence, not a real finding.")
    cols2 = ["substitution", "gene", "position",
             "BD_2017_2019_n", "BD_2017_2019_of",
             "BD_2021_2024_n", "BD_2021_2024_of",
             "BD_2025_n", "BD_2025_of"]
    print(chg[cols2].to_string(index=False))

    # ---- 3. where in the proteome ----
    print("\n== Substitutions by gene, Bangladeshi sequences ==")
    bd_any = tab[tab[[f"{c}_n" for c in bd_cols]].sum(axis=1) > 0]
    print(bd_any["gene"].value_counts().to_string())

    # ---- 4. envelope protein, annotated ----
    e = bd_any[bd_any["gene"].str.upper() == "E"].copy()
    if len(e):
        ann = e["position"].map(e_annotation)
        e["domain"] = [a for a, _ in ann]
        e["feature"] = [b for _, b in ann]
        e.to_csv(OUTDIR / "e_protein_domains.tsv", sep="\t", index=False)
        print(f"\n== Envelope protein substitutions ({len(e)}), annotated ==")
        print(e[["substitution", "position", "domain", "feature",
                 "BD_2017_2019_freq", "BD_2025_freq"]].to_string(index=False))
        print("\nBy domain:")
        print(e["domain"].value_counts().to_string())
        flagged = e[e["feature"] != ""]
        if len(flagged):
            print("\nSubstitutions at annotated functional sites "
                  "(worth checking against the antibody literature, not evidence of escape on their own):")
            print(flagged[["substitution", "domain", "feature"]].to_string(index=False))
    else:
        print("\nNo envelope substitutions found in the Bangladeshi sequences.")


if __name__ == "__main__":
    main()
