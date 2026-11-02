"""
08_curation_audit.py
--------------------
Audit the provenance of the data, which the pipeline so far has taken on trust.

Everything up to here has treated GenBank metadata as fact. This script asks the
awkward questions instead, because two of the project's claims currently rest on
unexamined provenance:

  Q1  Are the eighteen identical 2021 3III_A genomes one sequencing batch from one
      lab on one day? If so they are closer to ONE observation than eighteen, and
      "3III_A reached a third of 2021 samples" has to be reworded.

  Q2  Do the six 2025 whole genomes fail in the SAME places? All are exactly
      10,461 bases, which is a pipeline signature. If their missing regions
      coincide, a substitution shared by all six may be an assembly artefact
      rather than biology.

  Q3  Is any 2025 substitution supported by BOTH 2025 batches (the 10,461-base
      whole genomes AND the ~2,154-base amplicons)? Cross-batch support is much
      harder to explain as artefact, so this is the test that rescues or sinks
      the prM result.

  Q4  Are there probable duplicate patients or re-sequenced isolates - same
      country, date, division, submitter and length? We only de-duplicated
      identical sequences; near-identical ones from one patient still count twice.

  Q5  Are any dates impossible or suspicious - collected after release, in the
      future, or a whole batch sharing one date (a placeholder)?

  Q6  Why was MN018377.1 (China 2015, the earliest whole genome of the lineage)
      dropped? It is the sequence that would most help date the origin.

  Q7  How concentrated is each year's Bangladeshi data by submitter and division?
      Six 2025 genomes all from Dhaka is not the same as "Bangladesh".

Nothing is corrected here. The point is to see the problems and then decide,
on the record, what to do about each one.

Inputs : results/metadata/metadata.tsv
         results/metadata/metadata_clades.tsv
         results/nextclade/nextclade.tsv
Output : results/audit/*.tsv  plus a printed report

Run from the project folder:
  python scripts/08_curation_audit.py
"""

import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

MD = Path("results/metadata/metadata.tsv")
MDC = Path("results/metadata/metadata_clades.tsv")
NC = Path("results/nextclade/nextclade.tsv")
OUT = Path("results/audit")

INCURSION_2021 = "3III_A"
RESCUE = "MN018377.1"          # China 2015, dropped by the frameshift filter

# The two 2025 batches are defined by the data, not by accession prefixes.
# (An earlier version hardcoded prefixes and silently missed OZ484493.1,
# because the six accessions do not share a 6-character prefix.)
WG_MIN_LEN = 10000             # the 10,461-base whole-genome batch
BATCH_YEAR = 2025


def hr(title):
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72)


def parse_ranges(cell):
    """'996-1566,1942-2138,3658' -> [(996,1566),(1942,2138),(3658,3658)]"""
    out = []
    if not isinstance(cell, str) or not cell.strip() or cell.strip().lower() == "nan":
        return out
    for tok in cell.split(","):
        tok = tok.strip()
        m = re.match(r"^(\d+)-(\d+)$", tok)
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
            continue
        m = re.match(r"^(\d+)$", tok)
        if m:
            p = int(m.group(1))
            out.append((p, p))
    return out


def main():
    for p in (MD, MDC, NC):
        if not p.exists():
            sys.exit(f"Input not found: {p}")
    OUT.mkdir(parents=True, exist_ok=True)

    md = pd.read_csv(MD, sep="\t", low_memory=False)
    mdc = pd.read_csv(MDC, sep="\t", low_memory=False)
    nc = pd.read_csv(NC, sep="\t", low_memory=False)
    nc["accession"] = nc["seqName"].astype(str).str.split().str[0]

    d = mdc.merge(md[["accession", "submitter_affiliation", "release_date"]]
                  .drop_duplicates("accession"),
                  on="accession", how="left", suffixes=("", "_md"))
    for c in ("submitter_affiliation", "release_date"):
        if c not in d.columns and f"{c}_md" in d.columns:
            d[c] = d[f"{c}_md"]
    d["year"] = pd.to_numeric(d["year"], errors="coerce")
    bd = d[d["country"] == "Bangladesh"].copy()

    # ---------------------------------------------------------------- Q1
    hr("Q1. The eighteen identical 2021 3III_A genomes: one batch or eighteen?")
    inc = bd[(bd["clade"] == INCURSION_2021)].copy()
    print(f"  sequences: {len(inc)}")
    if len(inc):
        cols = [c for c in ["accession", "date", "division", "length",
                            "release_date", "submitter_affiliation"] if c in inc.columns]
        print(inc[cols].sort_values("accession").to_string(index=False))
        print("\n  submitters:")
        print("   " + inc["submitter_affiliation"].fillna("(none)").value_counts()
              .to_string().replace("\n", "\n   "))
        print("  release dates:")
        print("   " + inc["release_date"].fillna("(none)").value_counts()
              .to_string().replace("\n", "\n   "))
        print("  collection dates:")
        print("   " + inc["date"].fillna("(none)").value_counts()
              .to_string().replace("\n", "\n   "))
        n_sub = inc["submitter_affiliation"].nunique(dropna=True)
        n_rel = inc["release_date"].nunique(dropna=True)
        n_dat = inc["date"].nunique(dropna=True)
        print()
        if n_sub <= 1 and n_rel <= 1:
            print("  VERDICT: one submitter, one release date -> ONE BATCH.")
            print("  Eighteen identical genomes from one batch are close to ONE")
            print("  observation. Reword any claim about how common 3III_A was in 2021.")
        else:
            print(f"  VERDICT: {n_sub} submitters, {n_rel} release dates, "
                  f"{n_dat} collection dates.")
            print("  Identical genomes across more than one batch are much more likely")
            print("  to be a real tight transmission cluster.")
        inc.to_csv(OUT / "q1_2021_incursion.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- Q2
    hr("Q2. Do the six 2025 whole genomes fail in the same places?")
    miss_col = "missing" if "missing" in nc.columns else None
    b25 = bd[bd["year"] == BATCH_YEAR]
    wg25 = sorted(b25[b25["length"] >= WG_MIN_LEN]["accession"])
    amp25 = sorted(b25[b25["length"] < WG_MIN_LEN]["accession"])
    print(f"  whole-genome batch ({len(wg25)}): {', '.join(map(str, wg25))}")
    print(f"  amplicon batch ({len(amp25)}): {', '.join(map(str, amp25))}\n")
    if miss_col and wg25:
        per = {}
        for a in wg25:
            row = nc[nc["accession"] == a]
            per[a] = parse_ranges(row.iloc[0][miss_col]) if len(row) else []
        for a in wg25:
            tot = sum(b - x + 1 for x, b in per[a])
            print(f"  {a}  {tot:5d} missing bases in {len(per[a])} blocks")
            print(f"      {', '.join(f'{x}-{b}' for x, b in per[a][:12])}")
        # how often is each base missing across the six?
        cov = Counter()
        for a in wg25:
            for x, b in per[a]:
                for pos in range(x, b + 1):
                    cov[pos] += 1
        shared = sorted(p for p, n in cov.items() if n == len(wg25))
        print(f"\n  bases missing in ALL {len(wg25)} genomes: {len(shared)}")
        if shared:
            blocks, start, prev = [], shared[0], shared[0]
            for p in shared[1:]:
                if p != prev + 1:
                    blocks.append((start, prev)); start = p
                prev = p
            blocks.append((start, prev))
            print("  shared missing blocks: " +
                  ", ".join(f"{a}-{b}" for a, b in blocks[:20]))
            print("\n  VERDICT: the failures overlap. A substitution called in all six")
            print("  genomes OUTSIDE these blocks is fine; anything adjacent to them")
            print("  needs the cross-batch check in Q3 before it is believed.")
        else:
            print("  VERDICT: no base is missing in all six - the failures are not")
            print("  identical, which argues against a single systematic artefact.")
        pd.DataFrame({"position": sorted(cov), "n_genomes_missing":
                      [cov[p] for p in sorted(cov)]}).to_csv(
            OUT / "q2_2025_missing_overlap.tsv", sep="\t", index=False)
    else:
        print("  no 'missing' column, or no 2025 whole genomes found")

    # ---------------------------------------------------------------- Q3
    hr("Q3. Is any 2025 substitution supported by BOTH batches?")
    sub_col = "aaSubstitutions" if "aaSubstitutions" in nc.columns else None
    if sub_col:
        def subs_of(acc):
            r = nc[nc["accession"] == acc]
            if not len(r):
                return set()
            cell = r.iloc[0][sub_col]
            if not isinstance(cell, str):
                return set()
            return {t.strip() for t in cell.split(",") if ":" in t}

        wg_sets = {a: subs_of(a) for a in wg25}
        amp_sets = {a: subs_of(a) for a in amp25}
        print(f"  whole-genome batch: {len(wg25)}   amplicon batch: {len(amp25)}")

        in_all_wg = set.intersection(*wg_sets.values()) if wg_sets else set()
        in_all_amp = set.intersection(*amp_sets.values()) if amp_sets else set()
        both = sorted(in_all_wg & in_all_amp)
        only_wg = sorted(in_all_wg - in_all_amp)
        only_amp = sorted(in_all_amp - in_all_wg)

        print(f"\n  in every whole genome AND every amplicon ({len(both)}) "
              f"- CROSS-BATCH SUPPORTED, hardest to dismiss:")
        print("   " + (", ".join(both) if both else "(none)"))
        print(f"\n  in every whole genome only ({len(only_wg)}) - could be that "
              f"batch's pipeline, or a region the amplicons do not cover:")
        print("   " + (", ".join(only_wg) if only_wg else "(none)"))
        print(f"\n  in every amplicon only ({len(only_amp)}):")
        print("   " + (", ".join(only_amp) if only_amp else "(none)"))
        pd.DataFrame({"substitution": both + only_wg + only_amp,
                      "support": (["both batches"] * len(both)
                                  + ["whole-genome batch only"] * len(only_wg)
                                  + ["amplicon batch only"] * len(only_amp))}
                     ).to_csv(OUT / "q3_2025_cross_batch.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- Q4
    hr("Q4. Probable duplicate patients or re-sequenced isolates")
    key = [c for c in ["country", "date", "division", "length",
                       "submitter_affiliation"] if c in d.columns]
    grp = d.groupby(key, dropna=False)["accession"].apply(list)
    dups = {k: v for k, v in grp.items() if len(v) > 1}
    print(f"  groups sharing {' + '.join(key)}: {len(dups)}")
    bd_dups = {k: v for k, v in dups.items() if k[0] == "Bangladesh"}
    print(f"  of those, Bangladeshi: {len(bd_dups)}")
    shown = 0
    for k, v in sorted(bd_dups.items(), key=lambda kv: -len(kv[1])):
        print(f"    {len(v):3d}  {k}")
        print(f"         {', '.join(map(str, v[:10]))}"
              f"{' ...' if len(v) > 10 else ''}")
        shown += 1
        if shown >= 12:
            print("    (truncated)")
            break
    print("\n  These are NOT necessarily duplicates - a hospital can sequence many")
    print("  patients on one day. But each group is a candidate, and the ones that")
    print("  also have near-identical sequence deserve a look before they are")
    print("  counted as independent samples.")
    pd.DataFrame([{"group": str(k), "n": len(v), "accessions": ";".join(map(str, v))}
                  for k, v in dups.items()]).to_csv(
        OUT / "q4_possible_duplicates.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- Q5
    hr("Q5. Impossible or suspicious dates")
    dd = d.dropna(subset=["date"]).copy()
    dd["coll"] = pd.to_datetime(dd["date"], errors="coerce", format="mixed")
    dd["rel"] = pd.to_datetime(dd.get("release_date"), errors="coerce", format="mixed")
    after = dd[(dd["coll"].notna()) & (dd["rel"].notna()) & (dd["coll"] > dd["rel"])]
    print(f"  collected AFTER release: {len(after)}")
    if len(after):
        print(after[["accession", "country", "date", "release_date"]]
              .head(15).to_string(index=False))
    future = dd[dd["coll"] > pd.Timestamp.today()]
    print(f"  collection date in the future: {len(future)}")
    if len(future):
        print(future[["accession", "country", "date"]].head(15).to_string(index=False))

    print("\n  Bangladeshi collection dates shared by 3+ sequences "
          "(a shared date can be a batch placeholder):")
    vc = bd["date"].value_counts()
    rep = vc[vc >= 3]
    if len(rep):
        for dt, n in rep.items():
            accs = bd[bd["date"] == dt]["accession"].tolist()
            subs = bd[bd["date"] == dt]["submitter_affiliation"].nunique(dropna=True)
            flag = "  <-- one submitter: likely a batch date" if subs <= 1 else ""
            print(f"    {dt}: {n} sequences, {subs} submitter(s){flag}")
            print(f"       {', '.join(map(str, accs[:10]))}"
                  f"{' ...' if len(accs) > 10 else ''}")
    else:
        print("    none")

    # ---------------------------------------------------------------- Q6
    hr(f"Q6. Why was {RESCUE} dropped?")
    r = mdc[mdc["accession"] == RESCUE]
    if len(r):
        r = r.iloc[0]
        for c in ["country", "date", "length", "clade", "qc.overallStatus",
                  "totalMissing", "qc.frameShifts.status", "qc.stopCodons.status"]:
            if c in mdc.columns:
                print(f"  {c:26s} = {r[c]}")
        nr = nc[nc["accession"] == RESCUE]
        if len(nr):
            for c in nc.columns:
                if any(k in c.lower() for k in ("frameshift", "stopcodon", "coverage")):
                    print(f"  {c:26s} = {str(nr.iloc[0][c])[:200]}")
        print("\n  A frameshift in ONE gene does not make the rest of the genome")
        print("  unusable. If the frameshift is outside the region used for the tree,")
        print("  this sequence can be kept - and as the earliest whole genome of the")
        print("  lineage it is the single most useful sequence for dating the origin.")
    else:
        print(f"  {RESCUE} is not in metadata_clades.tsv at all.")

    # ---------------------------------------------------------------- Q7
    hr("Q7. How concentrated is the Bangladeshi data?")
    for col in ("division", "submitter_affiliation"):
        if col not in bd.columns:
            continue
        print(f"\n  by year x {col}:")
        t = pd.crosstab(bd["year"].fillna(0).astype(int), bd[col].fillna("(none)"))
        print("   " + t.to_string().replace("\n", "\n   "))
    print("\n  If a year's sequences come from one division and one submitter, that")
    print("  year describes that catchment, not the country. Say so in the paper.")

    print(f"\n\nTables written to {OUT}/")
    print("""
Decide and record, for each question: keep as is / reword the claim / rerun.
An audit that changes nothing is still worth having on the record - but an audit
whose findings are not acted on is worse than none.
""")


if __name__ == "__main__":
    main()
