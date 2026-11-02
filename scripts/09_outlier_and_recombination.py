"""
09_outlier_and_recombination.py
-------------------------------
Why is PQ815901.1 so divergent, and is anything in this dataset recombinant?

PQ815901.1 (Bangladesh, 2023) has a root-to-tip distance of 0.0664 while the 2025
Bangladeshi genomes sit at 0.0543-0.0552 and its sister OR821935.1 (USA 2023) sits
at 0.0521. Sister taxa share their whole path to the root, so roughly 0.014
substitutions per site - about 150 changes - are on its terminal branch alone.

Three explanations, and they leave different fingerprints:

  errors         extra differences scattered uniformly across the genome
  recombination  extra differences CONCENTRATED in one stretch, where the
                 sequence matches a DIFFERENT lineage than it does elsewhere
  fast evolution only credible once the first two are excluded

So this script does two things:

  1. For every sequence, count differences from the consensus of its own lineage
     and test whether they are clustered or spread (a dispersion test). A clustered
     excess is the recombination signature.

  2. A sliding-window best-match scan: in each window, which lineage consensus is
     this sequence closest to? A sequence that matches lineage A across most of the
     genome and lineage B in one block is a recombinant, and the window where the
     match switches is the approximate breakpoint.

This is a screen, not a substitute for GARD or 3SEQ. It is designed to tell you
quickly whether a formal recombination analysis is needed and where to look.

Inputs : results/datasets/wg.fasta
         results/datasets/wg_traits.tsv
Outputs: results/recombination/per_sequence.tsv      divergence + clustering per sequence
         results/recombination/window_match.tsv      best-matching lineage per window
         results/recombination/<acc>_profile.tsv     per-window detail for flagged ones

Run from the project folder:
  python scripts/09_outlier_and_recombination.py
  python scripts/09_outlier_and_recombination.py --focus PQ815901.1 --window 300
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

DATA = Path("results/datasets")
OUT = Path("results/recombination")
BASES = set("ACGT")


def load(stem):
    seqs = {}
    with (DATA / f"{stem}.fasta").open() as fh:
        for h, s in SimpleFastaParser(fh):
            seqs[h.split()[0]] = s.upper()
    traits = pd.read_csv(DATA / f"{stem}_traits.tsv", sep="\t")
    return seqs, traits


def consensus(seqs, names):
    """Majority base per column, ignoring gaps and Ns. '-' where nobody has data."""
    if not names:
        return ""
    L = len(seqs[names[0]])
    arr = [seqs[n] for n in names if n in seqs]
    out = []
    for i in range(L):
        c = Counter(s[i] for s in arr if s[i] in BASES)
        out.append(c.most_common(1)[0][0] if c else "-")
    return "".join(out)


def diffs(seq, ref):
    """Positions where both have a real base and they differ."""
    return [i for i, (a, b) in enumerate(zip(seq, ref))
            if a in BASES and b in BASES and a != b]


def covered_positions(seq):
    return sum(1 for c in seq if c in BASES)


def dispersion(positions, length, nwin):
    """
    Index of dispersion (variance / mean) of counts per window.
    ~1 means the differences are spread randomly (Poisson-like) - the error
    signature. Much greater than 1 means clustered - the recombination signature.
    """
    if not positions:
        return float("nan"), []
    edges = np.linspace(0, length, nwin + 1)
    counts, _ = np.histogram(positions, bins=edges)
    mean = counts.mean()
    if mean == 0:
        return float("nan"), counts.tolist()
    return float(counts.var() / mean), counts.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default="wg")
    ap.add_argument("--window", type=int, default=500,
                    help="sliding window size in bases for the match scan")
    ap.add_argument("--step", type=int, default=250)
    ap.add_argument("--focus", default="PQ815901.1",
                    help="sequence to profile in detail")
    ap.add_argument("--min-group", type=int, default=3,
                    help="minimum sequences for a lineage to get a consensus")
    args = ap.parse_args()

    if not (DATA / f"{args.stem}.fasta").exists():
        sys.exit(f"Missing {DATA}/{args.stem}.fasta - run scripts/04_select_sequences.py")
    seqs, traits = load(args.stem)
    OUT.mkdir(parents=True, exist_ok=True)
    L = len(next(iter(seqs.values())))
    print(f"{len(seqs)} sequences, alignment length {L}")

    # Lineage consensuses. Bangladeshi 3I_A.2 is split out because that is the
    # group the focus sequence belongs to, and a sequence must not be compared
    # against a consensus it helped define.
    traits = traits[traits["name"].isin(seqs)].copy()
    traits["group"] = np.where(
        (traits["clade"] == "3I_A.2") & (traits["country"] == "Bangladesh"),
        "BD_3I_A.2", traits["clade"])
    groups = {g: sorted(v) for g, v in traits.groupby("group")["name"]
              if len(v) >= args.min_group}
    print("\nGroups used for consensus:")
    for g, names in groups.items():
        print(f"  {g:12s} {len(names):4d}")

    cons = {g: consensus(seqs, names) for g, names in groups.items()}

    # ---------------- 1. per-sequence divergence and clustering -------------
    nwin_disp = 40
    rows = []
    for name, seq in seqs.items():
        g = traits.loc[traits["name"] == name, "group"]
        g = g.iloc[0] if len(g) else None
        if g not in cons:
            continue
        # rebuild the consensus without this sequence if the group is small,
        # so a sequence is never compared against itself
        others = [n for n in groups[g] if n != name]
        ref = cons[g] if len(others) >= args.min_group else consensus(seqs, others)
        if not ref:
            continue
        pos = diffs(seq, ref)
        covd = covered_positions(seq)
        disp, counts = dispersion(pos, L, nwin_disp)
        rows.append({
            "name": name, "group": g,
            "covered_bases": covd,
            "diffs_from_group_consensus": len(pos),
            "diffs_per_kb": round(1000 * len(pos) / covd, 3) if covd else None,
            "dispersion_index": round(disp, 2) if disp == disp else None,
            "max_window_count": max(counts) if counts else 0,
        })
    per = pd.DataFrame(rows).sort_values("diffs_per_kb", ascending=False)
    per.to_csv(OUT / "per_sequence.tsv", sep="\t", index=False)

    print("\n" + "=" * 72)
    print("  Most divergent sequences from their own lineage consensus")
    print("=" * 72)
    print(per.head(15).to_string(index=False))
    print("""
  dispersion_index ~1   differences spread evenly  -> looks like sequencing error
  dispersion_index >>1  differences clustered      -> looks like recombination""")

    # ---------------- 2. sliding-window best match --------------------------
    print("\n" + "=" * 72)
    print(f"  Window scan for {args.focus} (window {args.window}, step {args.step})")
    print("=" * 72)
    focus = args.focus
    if focus not in seqs:
        print(f"  {focus} is not in this dataset - skipping the window scan")
    else:
        fseq = seqs[focus]
        own = traits.loc[traits["name"] == focus, "group"]
        own = own.iloc[0] if len(own) else None
        # consensus of the focus sequence's own group, excluding itself
        own_ref = consensus(seqs, [n for n in groups.get(own, []) if n != focus])
        refs = {g: (own_ref if g == own else c) for g, c in cons.items()}

        prof = []
        for start in range(0, L - args.window + 1, args.step):
            end = start + args.window
            sub = fseq[start:end]
            if sum(1 for c in sub if c in BASES) < args.window * 0.5:
                continue
            d = {}
            for g, r in refs.items():
                if not r:
                    continue
                rs = r[start:end]
                comp = [(a, b) for a, b in zip(sub, rs) if a in BASES and b in BASES]
                if len(comp) < args.window * 0.5:
                    continue
                d[g] = sum(1 for a, b in comp if a != b) / len(comp)
            if not d:
                continue
            best = min(d, key=d.get)
            row = {"start": start + 1, "end": end, "best_match": best,
                   "best_dist": round(d[best], 4)}
            for g in refs:
                row[f"d_{g}"] = round(d[g], 4) if g in d else None
            prof.append(row)

        pdf = pd.DataFrame(prof)
        pdf.to_csv(OUT / f"{focus}_profile.tsv", sep="\t", index=False)
        if len(pdf):
            print(pdf.to_string(index=False))
            switches = pdf["best_match"].ne(pdf["best_match"].shift()).sum() - 1
            print(f"\n  best-match lineage changes {switches} time(s) across the genome")
            print("  " + " -> ".join(pdf["best_match"].loc[
                pdf["best_match"].ne(pdf["best_match"].shift())].tolist()))
            if switches == 0:
                print("""
  VERDICT: the same lineage is the closest match everywhere. No recombination
  signal. The excess divergence is therefore spread across the genome, which
  points to sequencing or assembly error rather than recombination - check the
  dispersion index above, and treat the sequence's private substitutions with
  caution.""")
            else:
                print("""
  VERDICT: the closest lineage CHANGES along the genome. That is the
  recombination signature. Note the windows where it switches - those bracket
  the breakpoint - and confirm with a formal method (GARD or 3SEQ) before
  reporting it. A recombinant must be excluded from the dated tree, because
  one tree cannot describe a genome with two histories.""")

    # ---------------- 3. whole-dataset switch scan --------------------------
    print("\n" + "=" * 72)
    print("  Any other sequence whose best match changes along the genome?")
    print("=" * 72)
    flagged = []
    for name, seq in seqs.items():
        g = traits.loc[traits["name"] == name, "group"]
        g = g.iloc[0] if len(g) else None
        if g not in cons:
            continue
        seq_best = []
        for start in range(0, L - args.window + 1, args.window):
            end = start + args.window
            sub = seq[start:end]
            if sum(1 for c in sub if c in BASES) < args.window * 0.5:
                continue
            d = {}
            for gg, r in cons.items():
                rs = r[start:end]
                comp = [(a, b) for a, b in zip(sub, rs) if a in BASES and b in BASES]
                if len(comp) < args.window * 0.5:
                    continue
                d[gg] = sum(1 for a, b in comp if a != b) / len(comp)
            if d:
                seq_best.append(min(d, key=d.get))
        if not seq_best:
            continue
        runs = [k for i, k in enumerate(seq_best) if i == 0 or k != seq_best[i - 1]]
        if len(runs) > 1:
            flagged.append({"name": name, "group": g, "n_switches": len(runs) - 1,
                            "pattern": " -> ".join(runs[:8])})
    fl = pd.DataFrame(flagged)
    fl.to_csv(OUT / "window_match.tsv", sep="\t", index=False)
    if len(fl):
        print(fl.sort_values("n_switches", ascending=False).to_string(index=False))
        print("\n  A sequence assigned to its own group everywhere is clean. Switching")
        print("  is a candidate recombinant - or simply a sequence near the boundary")
        print("  between two closely related lineage consensuses, which is common and")
        print("  harmless. Judge by how many windows switch and how consistent they are.")
    else:
        print("  None. Every sequence matches one lineage consistently across the")
        print("  genome, so there is no recombination signal in this dataset.")

    print(f"\nTables written to {OUT}/")


if __name__ == "__main__":
    main()
