"""
17_substitution_provenance.py
-----------------------------
For ONE substitution, show every Bangladeshi sequence individually: does it
carry the substitution, does it merely cover the position, or is it blind there.
Then cross-tabulate by year, by sequencing batch, and by submitting institution.

Why this script exists:

  results/audit/q3_2025_cross_batch.tsv cannot answer this. It was built with
  set.intersection, so it only lists substitutions present in EVERY sequence of
  a batch - substitutions that are FIXED. A substitution sitting at 36% can
  never appear in it, however real it is. Looking for a 36% substitution in that
  file returns nothing and that nothing means nothing.

  The question that actually matters for a substitution that ROSE over time is:
  are the carriers spread across years and across independent submissions, or
  are they one submission that happened to land late in the series? One
  submission is a laboratory. Several, across years, is a virus.

A substitution passes this check when:
  * carriers appear in more than one year, AND
  * carriers appear in more than one submitting institution or batch, AND
  * the non-carriers from the same years and batches are genuinely covered at
    that position (otherwise "absent early" just means "not sequenced early")

The third point is the one people forget. If every early sequence is blind at
the position, a 0/83 is not a zero.

Input : results/nextclade/nextclade.tsv
        results/metadata/metadata_clades.tsv
        results/nextclade/dataset/genome_annotation.gff3
Output: results/mutations/provenance_<SUB>.tsv  (one row per sequence)
        printed cross-tabs

Run from the project folder:
  python scripts/17_substitution_provenance.py --sub E:V139I
  python scripts/17_substitution_provenance.py --sub prM:N69K
  python scripts/17_substitution_provenance.py --sub E:V139I --country Bangladesh
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

NC = Path("results/nextclade/nextclade.tsv")
MD = Path("results/metadata/metadata_clades.tsv")
GFF = Path("results/nextclade/dataset/genome_annotation.gff3")
OUTDIR = Path("results/mutations")


def parse_unknown_ranges(cell):
    cov = defaultdict(list)
    if not isinstance(cell, str) or not cell.strip():
        return cov
    for tok in cell.split(","):
        tok = tok.strip()
        m = re.match(r"^([A-Za-z0-9_]+):(\d+)-(\d+)$", tok)
        if m:
            cov[m.group(1)].append((int(m.group(2)), int(m.group(3)))); continue
        m = re.match(r"^([A-Za-z0-9_]+):(\d+)$", tok)
        if m:
            p = int(m.group(2)); cov[m.group(1)].append((p, p))
    return cov


def read_gene_coords(gff_path):
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
    """Same test as 05_mutation_profile.py: in the aligned region AND not an N."""
    gc = gene_coords.get(gene)
    if gc is None or pd.isna(aln_start) or pd.isna(aln_end):
        return False
    codon_start = gc[0] + (pos - 1) * 3
    if codon_start < aln_start or codon_start + 2 > aln_end:
        return False
    for a, b in unknown.get(gene, []):
        if a <= pos <= b:
            return False
    return True


def batch_label(length):
    """Crude but honest: what kind of sequencing run produced this record."""
    if pd.isna(length):
        return "unknown length"
    n = int(length)
    if n >= 9000:
        return f"whole genome (~{round(n, -2)} nt)"
    if n >= 1200:
        return f"prM/E amplicon (~{round(n, -2)} nt)"
    return f"short fragment (~{round(n, -2)} nt)"


def crosstab(df, by, title):
    print(f"\n{title}")
    t = (df.groupby([by, "status"]).size().unstack(fill_value=0)
           .reindex(columns=["carrier", "covered, not carrier", "blind"],
                    fill_value=0))
    t["covered_total"] = t["carrier"] + t["covered, not carrier"]
    t["freq"] = [(f"{c / n:.2f}" if n else "-")
                 for c, n in zip(t["carrier"], t["covered_total"])]
    print(t.to_string())
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, help="e.g. E:V139I")
    ap.add_argument("--country", default="Bangladesh")
    ap.add_argument("--nextclade", default=str(NC))
    ap.add_argument("--metadata", default=str(MD))
    ap.add_argument("--gff", default=str(GFF))
    args = ap.parse_args()

    m = re.match(r"^([A-Za-z0-9_]+):([A-Z*\-])(\d+)([A-Z*\-])$", args.sub)
    if not m:
        sys.exit(f"--sub must look like E:V139I, got {args.sub!r}")
    gene, wt, pos, mut = m.group(1), m.group(2), int(m.group(3)), m.group(4)

    for p in (args.nextclade, args.metadata):
        if not Path(p).exists():
            sys.exit(f"missing input: {p}")
    nc = pd.read_csv(args.nextclade, sep="\t", low_memory=False)
    md = pd.read_csv(args.metadata, sep="\t", low_memory=False)

    gene_coords = read_gene_coords(Path(args.gff))
    if gene not in gene_coords:
        sys.exit(f"gene {gene!r} not in {args.gff}. Genes present: "
                 f"{sorted(gene_coords) or '(none - file missing or unparsed)'}")

    # nextclade's seqName carries the accession; match it to the metadata
    key = "accession" if "accession" in nc.columns else None
    if key is None:
        nc["accession"] = nc["seqName"].astype(str).str.split().str[0]
        key = "accession"
    cols = {c: f"{c}_nc" for c in ["alignmentStart", "alignmentEnd",
                                   "aaSubstitutions", "unknownAaRanges"]
            if c in nc.columns and c in md.columns}
    nc = nc.rename(columns=cols)
    aln_s = cols.get("alignmentStart", "alignmentStart")
    aln_e = cols.get("alignmentEnd", "alignmentEnd")
    subs_c = cols.get("aaSubstitutions", "aaSubstitutions")
    unk_c = cols.get("unknownAaRanges", "unknownAaRanges")

    keep = [c for c in ["accession", "country", "division", "date", "year",
                        "length", "submitter_affiliation", "clade"]
            if c in md.columns]
    d = md[keep].merge(nc[[key, aln_s, aln_e, subs_c, unk_c]],
                       on="accession", how="inner")
    if args.country:
        d = d[d["country"].fillna("") == args.country]
    if d.empty:
        sys.exit(f"no sequences for country={args.country!r}")

    rows = []
    for _, r in d.iterrows():
        unknown = parse_unknown_ranges(r.get(unk_c))
        cov = covered(gene, pos, unknown, r.get(aln_s), r.get(aln_e), gene_coords)
        cell = r.get(subs_c)
        carried = (isinstance(cell, str)
                   and args.sub in {t.strip() for t in cell.split(",")})
        if carried and not cov:
            # Nextclade called it, so the position was read: trust the call.
            cov = True
        status = "carrier" if carried else ("covered, not carrier" if cov
                                            else "blind")
        rows.append({
            "accession": r["accession"],
            "year": r.get("year"),
            "date": r.get("date"),
            "division": r.get("division"),
            "length": r.get("length"),
            "batch": batch_label(r.get("length")),
            "submitter": (str(r.get("submitter_affiliation"))[:55]
                          if pd.notna(r.get("submitter_affiliation")) else "(none)"),
            "clade": r.get("clade"),
            "status": status,
        })
    out = pd.DataFrame(rows).sort_values(["year", "accession"])

    OUTDIR.mkdir(parents=True, exist_ok=True)
    f = OUTDIR / f"provenance_{args.sub.replace(':', '_')}.tsv"
    out.to_csv(f, sep="\t", index=False)

    n_car = (out["status"] == "carrier").sum()
    n_cov = (out["status"] != "blind").sum()
    print(f"\n{'=' * 74}")
    print(f"  {args.sub}   ({gene} residue {pos}, {wt} -> {mut})   "
          f"country = {args.country}")
    print(f"{'=' * 74}")
    print(f"  sequences examined : {len(out)}")
    print(f"  covered at {gene}:{pos}      : {n_cov}")
    print(f"  carrying {args.sub:<12}: {n_car}"
          f"{f'  ({n_car / n_cov:.1%} of covered)' if n_cov else ''}")

    crosstab(out, "year", "By year  (freq = carriers / covered)")
    crosstab(out, "batch", "By sequencing batch")
    crosstab(out, "submitter", "By submitting institution")

    car = out[out["status"] == "carrier"]
    print(f"\nCarrier accessions ({len(car)}):")
    for _, r in car.iterrows():
        print(f"  {r['accession']:<14} {str(r['year']):<6} {r['batch']:<26} "
              f"{r['submitter']}")

    print(f"\nPer-sequence table -> {f}")

    # ---------------------------------------------------------------- verdict
    yrs = sorted({y for y in car["year"].dropna().tolist()})
    subm = sorted({s for s in car["submitter"].tolist() if s != "(none)"})
    batches = sorted(set(car["batch"].tolist()))
    print(f"\n{'-' * 74}\nWhat this supports\n{'-' * 74}")
    print(f"  carriers span {len(yrs)} year(s): "
          f"{', '.join(str(int(y)) for y in yrs) if yrs else '(none)'}")
    print(f"  carriers span {len(batches)} batch type(s): "
          f"{', '.join(batches) if batches else '(none)'}")
    print(f"  carriers span {len(subm)} named submitter(s): "
          f"{', '.join(subm) if subm else '(none named)'}")

    blind = out[out["status"] == "blind"]
    if len(blind):
        by = blind.groupby("year").size()
        print(f"\n  {len(blind)} sequence(s) are BLIND at this position, by year:")
        print("    " + ", ".join(f"{int(y)}: {n}" for y, n in by.items()
                                 if pd.notna(y)))
        print("    An 'absent' year is only meaningful if its sequences were")
        print("    covered. Check these before reading any 0/N as a true zero.")

    # Independence can come from either axis, and they are not equivalent.
    # Several years is the stronger signal: a fixed pipeline bias applies to
    # every run of that pipeline equally, so it cannot manufacture a frequency
    # that CHANGES over time. Several submitters in one year is weaker evidence
    # of a time trend but does rule out a single lab's assembly quirk.
    if len(yrs) >= 2 and (len(subm) >= 2 or len(batches) >= 2):
        v = ("strong. Carriers appear in several years AND several independent "
             "submissions or\n  batch types. Neither a single lab's pipeline nor "
             "a single run explains that.")
    elif len(yrs) >= 2:
        v = ("holds. Carriers appear in several years, and a fixed pipeline bias "
             "applies to\n  every run equally, so it cannot produce a frequency "
             "that changes over time.\n  Still worth naming the submitter in the "
             "paper, since it is one group's data.")
    elif len(subm) >= 2 or len(batches) >= 2:
        v = ("partial. Carriers come from independent submissions, so this is not "
             "one lab's\n  assembly quirk - but they sit in a single year, so this "
             "is not a time trend\n  either. It is a snapshot, and should be "
             "described as one.")
    else:
        v = ("weak. The carriers are one year from one submission. A laboratory "
             "or assembly\n  difference explains them as easily as the virus does. "
             "Do not put this on a\n  structure yet.")
    print(f"\n  VERDICT: {v}")
    print("\n  Provenance is necessary, not sufficient. It rules out the batch")
    print("  explanation. It does not show the substitution does anything.")


if __name__ == "__main__":
    main()
