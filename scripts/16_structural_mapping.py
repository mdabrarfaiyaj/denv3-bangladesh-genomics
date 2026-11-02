"""
16_structural_mapping.py
------------------------
Put the substitutions onto the protein structure and measure what can be measured.

For every substitution it reports:
  * whether the residue is present in the crystal at all (many are not: 1UZG is
    the E ECTODOMAIN, so the stem and transmembrane anchor, residues ~395-495,
    simply do not exist in it)
  * what residue the crystal actually has there. If the structure says GLU and
    the substitution says V139I, something is wrong - different numbering,
    different strain, or the wrong chain. This check is first because getting
    numbering wrong is the classic way structural mapping goes quietly wrong.
  * solvent accessibility, absolute and relative. A buried residue cannot be
    part of an antibody footprint.
  * distance to any ligand in the structure. For NS5 (4V0Q) the ligand is SAH,
    the product of the methyltransferase reaction, so this is the distance to
    the methyl-donor pocket.
  * every other residue within a cutoff, which is the structural neighbourhood
    the substitution actually perturbs.
  * domain, and for E whether it falls in the fusion loop or a glycosylation site.

What this does NOT do, deliberately:
  * predict stability changes (ddG). Predictors are unreliable on viral surface
    proteins and a number with no error bar invites overclaiming.
  * claim immune escape. A substitution in an epitope is a hypothesis. Escape is
    a neutralisation result, and there is no neutralisation data here.

Structures (download first; both are DENV-3, so no homology modelling):
  wget -P data/structures https://files.rcsb.org/download/1UZG.pdb   # E ectodomain
  wget -P data/structures https://files.rcsb.org/download/4V0Q.pdb   # full-length NS5 + SAH

Run from the project folder:
  python scripts/16_structural_mapping.py --protein E
  python scripts/16_structural_mapping.py --protein NS5
  python scripts/16_structural_mapping.py --protein E --subs E:V139I,E:T303A
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser, Selection
from Bio.PDB.SASA import ShrakeRupley

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# Theoretical maximum solvent accessibility, Tien et al. 2013 (PLoS ONE 8:e80635),
# used to turn absolute SASA into a relative, comparable figure.
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# Envelope protein, standard flavivirus E numbering. DI is discontinuous.
E_DOMAINS = [(1, 52, "domain I"), (53, 131, "domain II"), (132, 193, "domain I"),
             (194, 279, "domain II"), (280, 296, "domain I"),
             (297, 394, "domain III"), (395, 495, "stem / anchor")]
E_FEATURES = [(98, 110, "fusion loop"), (67, 67, "N-glycosylation site N67"),
              (153, 153, "N-glycosylation site N153"), (62, 80, "glycan loop")]

# NS5: methyltransferase then RNA-dependent RNA polymerase.
NS5_DOMAINS = [(1, 265, "methyltransferase (MTase)"),
               (266, 320, "MTase-RdRp linker"),
               (321, 900, "RNA-dependent RNA polymerase (RdRp)")]

DEFAULTS = {
    "E":   {"pdb": "data/structures/1UZG.pdb", "chain": "A"},
    "NS5": {"pdb": "data/structures/4V0Q.pdb", "chain": "A"},
}


def annotate(protein, pos):
    table = E_DOMAINS if protein == "E" else NS5_DOMAINS
    dom = next((n for a, b, n in table if a <= pos <= b), "outside annotated range")
    feats = ""
    if protein == "E":
        feats = "; ".join(n for a, b, n in E_FEATURES if a <= pos <= b)
    return dom, feats


def load_substitutions(protein, explicit):
    """Either an explicit list, or everything significant plus the lineage markers."""
    if explicit:
        return [s.strip() for s in explicit.split(",") if s.strip()]
    subs = []
    t = Path("results/mutations/temporal_test.tsv")
    if t.exists():
        d = pd.read_csv(t, sep="\t")
        d = d[(d["p_adj"] < 0.05) & (d["direction"] == "rose")]
        subs += d["substitution"].tolist()
        print(f"  {len(d)} substitution(s) that rose significantly over time")
    l = Path("results/mutations/lineage_defining.tsv")
    if l.exists():
        d = pd.read_csv(l, sep="\t")
        subs += d["substitution"].tolist()
        print(f"  {len(d)} lineage-defining substitution(s)")
    q = Path("results/audit/q3_2025_cross_batch.tsv")
    if q.exists():
        d = pd.read_csv(q, sep="\t")
        d = d[d["support"] == "both batches"]
        subs += d["substitution"].tolist()
        print(f"  {len(d)} cross-batch-supported 2025 substitution(s)")
    out, seen = [], set()
    for s in subs:
        if s.startswith(protein + ":") and s not in seen:
            seen.add(s); out.append(s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", choices=["E", "NS5"], required=True)
    ap.add_argument("--pdb", default=None)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--subs", default=None,
                    help="comma-separated, e.g. E:V139I,E:T303A. Default: the "
                         "significant and lineage-defining ones from results/")
    ap.add_argument("--contact-cutoff", type=float, default=8.0)
    args = ap.parse_args()

    pdb = Path(args.pdb or DEFAULTS[args.protein]["pdb"])
    chain_id = args.chain or DEFAULTS[args.protein]["chain"]
    if not pdb.exists():
        sys.exit(f"Structure not found: {pdb}\n"
                 f"Download it first:\n"
                 f"  mkdir -p data/structures\n"
                 f"  wget -P data/structures "
                 f"https://files.rcsb.org/download/"
                 f"{'1UZG' if args.protein == 'E' else '4V0Q'}.pdb")

    print(f"\n== {args.protein}: {pdb} chain {chain_id} ==")
    subs = load_substitutions(args.protein, args.subs)
    if not subs:
        sys.exit(f"\nNo {args.protein} substitutions to map. Pass --subs, or run "
                 f"scripts/15_temporal_substitutions.py and 05_mutation_profile.py.")
    print(f"\nMapping {len(subs)}: {', '.join(subs)}")

    structure = PDBParser(QUIET=True).get_structure("s", str(pdb))
    model = structure[0]
    if chain_id not in [c.id for c in model]:
        sys.exit(f"chain {chain_id} not in structure; chains present: "
                 f"{[c.id for c in model]}")
    chain = model[chain_id]

    print("\nComputing solvent accessibility (Shrake-Rupley) ...")
    ShrakeRupley().compute(model, level="R")

    # residues by sequence number, and any ligands (HETATM, excluding water)
    residues, ligands = {}, []
    for res in chain:
        het, num, _ = res.id
        if het == " ":
            residues[num] = res
        elif res.get_resname() != "HOH":
            ligands.append(res)
    for c in model:
        for res in c:
            if res.id[0] != " " and res.get_resname() != "HOH" and res not in ligands:
                ligands.append(res)
    span = (min(residues), max(residues)) if residues else (None, None)
    print(f"  {len(residues)} residues modelled, numbering {span[0]}-{span[1]}")
    if ligands:
        print(f"  ligands present: "
              f"{', '.join(sorted({l.get_resname() for l in ligands}))}")

    all_atoms = Selection.unfold_entities(model, "A")
    rows = []
    for sub in subs:
        _, _, change = sub.partition(":")
        wt, pos, mut = change[0], int(change[1:-1]), change[-1]
        dom, feat = annotate(args.protein, pos)
        r = {"substitution": sub, "wt": wt, "position": pos, "mut": mut,
             "domain": dom, "feature": feat}

        res = residues.get(pos)
        if res is None:
            r["in_structure"] = "NO - not modelled in this crystal"
            rows.append(r); continue

        actual = THREE_TO_ONE.get(res.get_resname(), "?")
        r["in_structure"] = "yes"
        r["residue_in_pdb"] = actual
        r["numbering_check"] = ("match" if actual == wt
                                else f"MISMATCH (pdb has {actual}, expected {wt})")

        sasa = getattr(res, "sasa", None)
        r["sasa_A2"] = round(sasa, 1) if sasa is not None else None
        if sasa is not None and actual in MAX_ASA:
            rel = 100 * sasa / MAX_ASA[actual]
            r["rel_sasa_pct"] = round(rel, 1)
            r["exposure"] = ("buried" if rel < 10 else
                             "partly buried" if rel < 35 else "exposed")

        if ligands:
            best, bestd = None, 1e9
            for lig in ligands:
                for a1 in res:
                    for a2 in lig:
                        dd = a1 - a2
                        if dd < bestd:
                            bestd, best = dd, lig.get_resname()
            r["nearest_ligand"] = best
            r["dist_to_ligand_A"] = round(bestd, 2)

        near = set()
        for a1 in res:
            for a2 in all_atoms:
                p2 = a2.get_parent()
                if p2 is res or p2.id[0] != " ":
                    continue
                if (a1 - a2) <= args.contact_cutoff:
                    near.add((p2.id[1], THREE_TO_ONE.get(p2.get_resname(), "?")))
        r["n_contacts"] = len(near)
        r["contacts"] = ";".join(f"{aa}{n}" for n, aa in sorted(near))
        rows.append(r)

    out = pd.DataFrame(rows)
    if "n_contacts" in out.columns:
        out["n_contacts"] = out["n_contacts"].astype("Int64")
    dest = Path("results/structure"); dest.mkdir(parents=True, exist_ok=True)
    f = dest / f"{args.protein}_mapping.tsv"
    out.to_csv(f, sep="\t", index=False)

    show = [c for c in ["substitution", "in_structure", "numbering_check", "domain",
                        "feature", "rel_sasa_pct", "exposure",
                        "nearest_ligand", "dist_to_ligand_A", "n_contacts"]
            if c in out.columns]
    print(f"\n{'=' * 78}")
    print(out[show].to_string(index=False))
    print(f"{'=' * 78}")
    print(f"\nFull table (including contact residues) -> {f}")

    # a PyMOL session so the result can be looked at, not just read
    pml = dest / f"{args.protein}_view.pml"
    sel = "+".join(str(r["position"]) for r in rows
                   if r.get("in_structure") == "yes")
    with pml.open("w") as fh:
        fh.write(f"load {pdb.resolve()}\nhide everything\nshow cartoon\n"
                 f"color grey80\nbg_color white\n")
        if sel:
            fh.write(f"select subs, chain {chain_id} and resi {sel}\n"
                     f"show spheres, subs\ncolor red, subs\n")
        if ligands:
            fh.write("show sticks, hetatm and not resn HOH\n"
                     "color yellow, hetatm and not resn HOH\n")
        fh.write("orient\n")
    print(f"PyMOL session        -> {pml}   (pymol {pml})")

    print("""
Reading this:

  * numbering_check is the first thing to look at. A MISMATCH means the
    structure's numbering does not line up with the substitution's, and every
    other number in that row is meaningless until it is resolved.

  * rel_sasa_pct under ~10% means the residue is buried. A buried substitution
    cannot be part of an antibody footprint, whatever else it might do.

  * dist_to_ligand_A for NS5 is the distance to SAH, which occupies the
    methyl-donor pocket. Within ~8 A is the active site neighbourhood; 30 A
    away is a different part of the protein.

  * 'not modelled in this crystal' is an honest answer, not a failure. 1UZG is
    the ectodomain only, so the stem and anchor are absent.

  * None of this demonstrates immune escape. The defensible sentence is
    "residue X is surface-exposed in domain Y, within Z A of <feature>" - and
    then a citation to whether anyone has tested that position experimentally.
""")


if __name__ == "__main__":
    main()
