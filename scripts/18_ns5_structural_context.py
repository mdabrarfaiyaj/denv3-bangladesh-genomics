#!/usr/bin/env python3
"""
18_ns5_structural_context.py   (v2)
-----------------------------------
NS5 structural context for the lineage-defining substitutions of 3I_A.2.

Why this exists separately from 16_structural_mapping.py
--------------------------------------------------------
Script 16 reports the distance to the NEAREST het group, which for H52R was
3.28 A to ACT. ACT is acetate from the crystallisation buffer. The distance to
SAH - the methyl donor, and the ligand that marks the methyltransferase active
site - was never printed. This script classifies het groups and reports the
distance to EACH meaningful one, never the nearest anything.

ZN and MG are NOT buffer. Flavivirus NS5 carries structural zinc and catalytic
magnesium; treating them as contaminants would be the opposite error.

v2 changes
----------
  * MTase catalytic tetrad corrected to K61-D146-K180-E216, VERIFIED against
    Zhao et al. 2015 PNAS 112:14834 (full-length DENV NS5 with cap-0 RNA).
    v1 used 181/217, which are the West Nile virus numbers shifted; the crystal
    printed V181 and M217 and caught it. Every tetrad distance from v1 is void.
  * SAM-binding-pocket residue list REMOVED. The SAH ligand is the pocket, and
    its coordinates are measured, not guessed.
  * Cap-binding-site residue list REMOVED. In 5DTO the M7G cap is present as a
    ligand, so the distance is measured directly.
  * M7G added as a meaningful het group (the 7-methyl-GTP cap).
  * Nucleic-acid chains are now measured. 5DTO has cap-0 viral RNA in chain B;
    distance from a substitution to the RNA is a real number and was the thing
    4V0Q could not give.
  * Features that remain unverified are PREFIXED 'UNVERIFIED' in the output so
    a number cannot be quoted without seeing the warning.

What this does NOT do
---------------------
  * predict ddG or stability change
  * claim an effect on capping, replication or fitness. Proximity to an active
    site is a hypothesis. An effect is an enzyme assay, and there is none here.

Run:
  python scripts/18_ns5_structural_context.py --pdb data/structures/5DTO.pdb
  python scripts/18_ns5_structural_context.py --pdb data/structures/4V0Q.pdb
  python scripts/18_ns5_structural_context.py --subs NS5:H52R
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

MAX_ASA = {  # Tien et al. 2013, PLoS ONE 8:e80635
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# --- nucleic acid ------------------------------------------------------------
NUCLEOTIDES = {"A", "C", "G", "U", "I",
               "DA", "DC", "DG", "DT", "DU",
               "M7G", "7MG", "GTP"}

# --- het group classification -------------------------------------------------
MEANINGFUL_HET = {
    "SAH": "S-adenosyl-L-homocysteine (MTase product; marks the SAM pocket)",
    "SAM": "S-adenosyl-L-methionine (methyl donor)",
    "M7G": "7-methyl-GTP (the RNA cap; marks the cap-binding site)",
    "7MG": "7-methyl-guanosine (cap)",
    "GTP": "GTP (cap analogue)",
    "ZN":  "structural zinc (NS5 RdRp zinc site)",
    "MG":  "magnesium (catalytic metal)",
    "MN":  "manganese (catalytic metal)",
}
BUFFER_HET = {
    "ACT": "acetate (buffer)", "GOL": "glycerol (cryoprotectant)",
    "SO4": "sulfate (buffer)", "PO4": "phosphate (buffer)",
    "EDO": "ethylene glycol (cryoprotectant)", "PEG": "polyethylene glycol",
    "PG4": "PEG fragment", "1PE": "PEG fragment", "MPD": "MPD (cryoprotectant)",
    "CL": "chloride ion", "NA": "sodium ion", "K": "potassium ion",
    "TRS": "Tris (buffer)", "DMS": "DMSO", "IMD": "imidazole", "CIT": "citrate",
    "FMT": "formate", "NO3": "nitrate", "IOD": "iodide", "BR": "bromide",
    "CA": "calcium (often adventitious)",
}

# Domain and subdomain boundaries.
# MTase 1-265, linker 266-272, RdRp 273-900 (Lim et al. 2016; Zhao et al. 2015).
# RdRp subdomains for DENV3: fingers 273-315, 416-496, 543-600; the rest of the
# N-terminal extension and the NLS region sit between the fingers segments
# (Yap et al. 2007 J Virol 81:4753, DENV3 RdRp at 1.85 A). Ordered most
# specific first; the first match wins.
NS5_DOMAINS = [
    (1, 265, "MTase"),
    (266, 272, "linker"),
    (273, 315, "RdRp fingers"),
    (316, 415, "RdRp NLS region"),
    (416, 496, "RdRp fingers"),
    (497, 542, "RdRp palm"),
    (543, 600, "RdRp fingers"),
    (601, 713, "RdRp palm"),
    (714, 781, "RdRp thumb"),
    (782, 809, "RdRp priming loop"),
    (810, 900, "RdRp thumb"),
]
# Coarse domain, used for the interdomain-contact count.
NS5_COARSE = [(1, 265, "MTase"), (266, 272, "linker"), (273, 900, "RdRp")]

# --- named functional features ------------------------------------------------
# Every entry below is VERIFIED against a named paper except where the key
# begins APPROX. v1 of this script used an unverified tetrad (181/217) and the
# crystal caught it; nothing goes in here now without a citation.
FEATURES = {
    # Zhao et al. 2015, PNAS 112:14834 (DENV full-length NS5 with cap-0 RNA).
    # The flavivirus literature carries several numberings for this motif
    # (WNV is K61-D146-K182-E218). These are the dengue numbers.
    "MTase catalytic tetrad (K-D-K-E)": [61, 146, 180, 216],
    # Yap et al. 2007, J Virol 81:4753 (DENV3 RdRp, 1.85 A).
    "RdRp motif C (GDD)": [662, 663, 664],
    "RdRp motif A (D533)": [533],
    # Lim et al. 2016, PLoS Pathog 12:e1005737; Yap et al. 2007.
    # W795/S796/I797 are core members of this loop.
    "priming loop (782-809)": list(range(782, 810)),
    # Lim et al. 2016: the alpha/beta NLS lines the NITD-434 pocket together
    # with motif B and the priming loop.
    "alpha/beta NLS (369-415)": list(range(369, 416)),
    "RdRp motif B (597-613)": list(range(597, 614)),
    "RdRp motif F (468-473)": list(range(468, 474)),
    # H18-H19 forms one wall of the de novo initiation chamber with the priming
    # loop. The 395-426 boundary comes from a ZIKV-numbered analysis, so it is
    # APPROXIMATE in DENV - flagged in the key so a distance cannot be quoted
    # without the caveat. The chamber exists in the pre-initiation
    # conformation only; the priming loop retracts during elongation.
    "APPROX H18-H19 chamber (395-426)": list(range(395, 427)),
}

DEFAULT_PDB = "data/structures/5DTO.pdb"


def domain_of(pos):
    """Fine-grained subdomain."""
    return next((n for a, b, n in NS5_DOMAINS if a <= pos <= b), "outside range")


def coarse_domain_of(pos):
    """MTase / linker / RdRp, for the interdomain-contact count."""
    return next((n for a, b, n in NS5_COARSE if a <= pos <= b), "outside range")


def min_dist(res_a, res_b):
    best = None
    for x in res_a:
        if x.element == "H":
            continue
        for y in res_b:
            if y.element == "H":
                continue
            d = x - y
            if best is None or d < best:
                best = d
    return best


def load_subs(explicit):
    if explicit:
        return [s.strip() for s in explicit.split(",") if s.strip()]
    p = Path("results/mutations/lineage_defining.tsv")
    if not p.exists():
        sys.exit(f"{p} not found; pass --subs NS5:H52R,...")
    d = pd.read_csv(p, sep="\t")
    subs = [s for s in d["substitution"] if str(s).startswith("NS5:")]
    if not subs:
        sys.exit("No NS5 substitutions in lineage_defining.tsv")
    return subs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default=DEFAULT_PDB)
    ap.add_argument("--chain", default="A", help="protein chain to map onto")
    ap.add_argument("--subs", default=None)
    ap.add_argument("--contact-cutoff", type=float, default=8.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    pdb = Path(args.pdb)
    if not pdb.exists():
        sys.exit(f"Structure not found: {pdb}\n"
                 f"  wget -P data/structures https://files.rcsb.org/download/5DTO.pdb")
    out = Path(args.out or
               f"results/structure/ns5_context_{pdb.stem}.tsv")

    subs = load_subs(args.subs)
    print(f"\n== NS5 structural context: {pdb} chain {args.chain} ==")
    print(f"Substitutions: {', '.join(subs)}")

    model = PDBParser(QUIET=True).get_structure("s", str(pdb))[0]
    chains = [c.id for c in model]
    if args.chain not in chains:
        sys.exit(f"chain {args.chain} not present; chains: {chains}")
    print(f"Chains in structure: {', '.join(chains)}")
    chain = model[args.chain]

    print("Computing solvent accessibility (Shrake-Rupley) ...")
    ShrakeRupley().compute(model, level="R")

    residues = {r.id[1]: r for r in chain
                if r.id[0] == " " and r.get_resname().strip() in THREE_TO_ONE}
    if residues:
        print(f"  {len(residues)} protein residues in chain {args.chain}, "
              f"numbering {min(residues)}-{max(residues)}")

    # --- sort everything that is not the protein chain -------------------------
    meaningful, buffer_hets, unknown, nucleic = [], [], [], []
    for c in model:
        for r in c:
            name = r.get_resname().strip()
            if name == "HOH":
                continue
            if c.id == args.chain and r.id[0] == " " and name in THREE_TO_ONE:
                continue
            entry = (name, c.id, r.id[1], r)
            if name in NUCLEOTIDES and r.id[0] == " ":
                nucleic.append(entry)
            elif name in MEANINGFUL_HET:
                meaningful.append(entry)
                if name in ("M7G", "7MG"):
                    nucleic.append(entry)   # the cap is part of the RNA too
            elif name in BUFFER_HET:
                buffer_hets.append(entry)
            elif r.id[0] != " ":
                unknown.append(entry)

    print("\nHet groups and chains:")
    print(f"  meaningful ligands : "
          f"{sorted({n for n, _, _, _ in meaningful}) or 'none'}")
    print(f"  buffer (ignored)   : "
          f"{sorted({n for n, _, _, _ in buffer_hets}) or 'none'}")
    if unknown:
        print(f"  UNCLASSIFIED - check these yourself: "
              f"{sorted({n for n, _, _, _ in unknown})}")
    if nucleic:
        by_chain = {}
        for n, cid, num, _ in nucleic:
            by_chain.setdefault(cid, []).append(n)
        for cid, names in by_chain.items():
            print(f"  nucleic acid chain {cid}: {len(names)} residues "
                  f"({''.join(x if len(x) == 1 else f'[{x}]' for x in names)})")
    else:
        print("  no nucleic acid in this structure "
              "(4V0Q has none; use 5DTO for RNA distances)")

    # --- feature numbering check ----------------------------------------------
    print("\nFeature numbering check:")
    for fname, positions in FEATURES.items():
        got = []
        for p in positions:
            r = residues.get(p)
            got.append(f"{THREE_TO_ONE.get(r.get_resname(), '?')}{p}" if r
                       else f"-{p}")
        flag = "  <-- NOT VERIFIED" if fname.startswith("UNVERIFIED") else ""
        print(f"  {fname:34s} {' '.join(got)}{flag}")
    print("  (tetrad should read K61 D146 K180 E216; motif C should read "
          "G662 D663 D664)")

    rows = []
    for sub in subs:
        _, _, change = sub.partition(":")
        wt, pos, mut = change[0], int(change[1:-1]), change[-1]
        r = {"substitution": sub, "wt": wt, "position": pos, "mut": mut,
             "domain": domain_of(pos), "structure": pdb.stem}

        res = residues.get(pos)
        if res is None:
            r["in_structure"] = "NO - not modelled"
            rows.append(r)
            continue

        actual = THREE_TO_ONE.get(res.get_resname(), "?")
        r["in_structure"] = "yes"
        r["residue_in_pdb"] = actual
        r["numbering_check"] = ("match" if actual == wt
                                else f"MISMATCH (pdb {actual}, expected {wt})")

        sasa = getattr(res, "sasa", None)
        if sasa is not None:
            r["sasa_A2"] = round(sasa, 1)
            if actual in MAX_ASA:
                rel = 100 * sasa / MAX_ASA[actual]
                r["rel_sasa_pct"] = round(rel, 1)
                r["exposure"] = ("buried" if rel < 10 else
                                 "partly buried" if rel < 35 else "exposed")

        for name, cid, resnum, lig in meaningful:
            d = min_dist(res, lig)
            if d is not None:
                r[f"dist_{name}_{cid}{resnum}_A"] = round(d, 2)

        # distance to the nucleic acid as a whole, and to its nearest residue
        if nucleic:
            best = min(((min_dist(res, nuc), n, cid, num)
                        for n, cid, num, nuc in nucleic),
                       key=lambda t: (t[0] is None, t[0]))
            if best[0] is not None:
                r["dist_RNA_A"] = round(best[0], 2)
                r["nearest_RNA_residue"] = f"{best[1]}{best[3]} (chain {best[2]})"
                r["RNA_contact"] = "yes" if best[0] <= 4.5 else "no"

        if buffer_hets:
            nb = min(((min_dist(res, l), n) for n, c, i, l in buffer_hets),
                     key=lambda t: (t[0] is None, t[0]))
            if nb[0] is not None:
                r["nearest_BUFFER_het"] = f"{nb[1]} ({nb[0]:.2f} A) - IGNORE"

        for fname, positions in FEATURES.items():
            if pos in positions:
                # The residue is a member of this feature. Measuring it against
                # itself gives 0.00 A, which reads like a distance and is not
                # one. Say so instead.
                r[f"dist_{fname}_A"] = "INSIDE"
                continue
            ds = [min_dist(res, residues[p]) for p in positions
                  if p in residues and p != pos]
            ds = [d for d in ds if d is not None]
            if ds:
                r[f"dist_{fname}_A"] = round(min(ds), 2)
            else:
                # every position of this feature is disordered in this crystal
                r[f"dist_{fname}_A"] = "not modelled"

        contacts = []
        for num, other in residues.items():
            if num == pos:
                continue
            d = min_dist(res, other)
            if d is not None and d <= args.contact_cutoff:
                contacts.append((num, THREE_TO_ONE.get(other.get_resname(), "?"),
                                 coarse_domain_of(num), domain_of(num)))
        contacts.sort()
        r["n_contacts"] = len(contacts)
        by_dom, by_sub = {}, {}
        for num, aa, dom, sub in contacts:
            by_dom.setdefault(dom, []).append(f"{aa}{num}")
            by_sub.setdefault(sub, []).append(f"{aa}{num}")
        for dom in ("MTase", "linker", "RdRp"):
            r[f"contacts_{dom}_n"] = len(by_dom.get(dom, []))
        r["contacts_by_subdomain"] = "; ".join(
            f"{k}: {len(v)}" for k, v in sorted(by_sub.items()))
        r["contacts"] = ";".join(f"{aa}{num}" for num, aa, _, _ in contacts)
        own = coarse_domain_of(pos)
        r["cross_domain_contacts"] = sum(len(v) for d, v in by_dom.items()
                                         if d != own)
        r["interface_residue"] = "yes" if r["cross_domain_contacts"] >= 3 else "no"
        rows.append(r)

    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)

    print("\n--- summary ---")
    for _, row in df.iterrows():
        print(f"\n{row['substitution']}  [{row.get('domain')}]")
        if row.get("in_structure") != "yes":
            print("  not modelled in this crystal")
            continue
        print(f"  numbering      : {row.get('numbering_check')}")
        print(f"  exposure       : {row.get('exposure')} "
              f"({row.get('rel_sasa_pct')}% relative)")
        if pd.notna(row.get("dist_RNA_A")):
            print(f"  RNA            : {row['dist_RNA_A']:.2f} A to "
                  f"{row['nearest_RNA_residue']}  "
                  f"[direct contact: {row['RNA_contact']}]")
        for k in row.index:
            if k.startswith("dist_") and k != "dist_RNA_A" and pd.notna(row[k]):
                v = row[k]
                if isinstance(v, str):
                    print(f"  {k[5:-2]:32s} {v}")
                else:
                    print(f"  {k[5:-2]:32s} {v:.2f} A")
        if pd.notna(row.get("nearest_BUFFER_het")):
            print(f"  {row['nearest_BUFFER_het']}")
        print(f"  contacts       : {row.get('n_contacts')} total "
              f"(MTase {row.get('contacts_MTase_n')}, "
              f"linker {row.get('contacts_linker_n')}, "
              f"RdRp {row.get('contacts_RdRp_n')})")
        print(f"    by subdomain : {row.get('contacts_by_subdomain')}")
        print(f"  interdomain    : {row.get('cross_domain_contacts')} outside "
              f"own domain -> interface residue: {row.get('interface_residue')}")

    print(f"\nWrote {out}")
    print("\nReminder: proximity is not function. Nothing here demonstrates an "
          "effect on capping, replication or fitness.")


if __name__ == "__main__":
    main()
