#!/usr/bin/env python3
"""
20_ns5_figure.py
----------------
Figures for the NS5 structural result.

Produces three things:

  1. docs/ns5_3d.html   - a self-contained interactive 3D viewer (3Dmol.js).
                          The PDB coordinates are embedded in the file, so it
                          works offline and needs no RCSB access. Host it on
                          GitHub Pages and link it from the README.
  2. results/figures/ns5_structure.pml - a PyMOL script for a publication PNG.
  3. results/figures/ns5_conservation.png - the conservation panel, matplotlib.

Why both static and interactive
-------------------------------
GitHub strips JavaScript from rendered README.md, so an interactive viewer
cannot live in the README itself. The pattern that works is: a static PNG in
the README, which is a clickable link to the interactive page on GitHub Pages.

Run:
  python scripts/20_ns5_figure.py
  python scripts/20_ns5_figure.py --pdb data/structures/4V0Q.pdb
"""

import argparse
import json
from pathlib import Path

import pandas as pd

# Residues to highlight, with the colour used in every panel.
HIGHLIGHT = {
    52:  {"sub": "H52R",  "color": "#e6194b",
          "note": "MTase-RdRp interface; 12.7 A from the K-D-K-E tetrad"},
    139: {"sub": "K139R", "color": "#3cb44b",
          "note": "surface, 16.7 A from SAH; no structural context"},
    404: {"sub": "T404A", "color": "#4363d8",
          "note": "de novo initiation chamber; 11.7 A from the priming loop"},
}

# Verified feature ranges (see scripts/18 for citations).
FEATURE_SELECTIONS = {
    "MTase":        (1, 265,   "#d9d9d9"),
    "linker":       (266, 320, "#bfbfbf"),
    "RdRp":         (321, 900, "#a6a6a6"),
}
MOTIF_C = (662, 664)
PRIMING_LOOP = (782, 809)
TETRAD = [61, 146, 180, 216]


def build_3dmol_html(pdb_path, out_path, title):
    """Self-contained 3Dmol.js page with the coordinates embedded."""
    pdb_text = Path(pdb_path).read_text()
    # strip the bulkiest records that add nothing to a cartoon view
    keep = [ln for ln in pdb_text.splitlines()
            if ln.startswith(("ATOM", "HETATM", "TER", "END", "HEADER",
                              "TITLE", "SEQRES", "HELIX", "SHEET"))]
    pdb_js = json.dumps("\n".join(keep))

    highlight_js = json.dumps(
        [{"resi": k, "color": v["color"], "label": v["sub"], "note": v["note"]}
         for k, v in HIGHLIGHT.items()])

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1a1a; --muted: #5a5a5a; --line: #e2e2e2;
    --panel: #fafafa;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --bg: #15171a; --fg: #e8e8e8; --muted: #a0a0a0; --line: #2c2f33;
      --panel: #1c1f23;
    }}
  }}
  :root[data-theme="dark"] {{
    --bg: #15171a; --fg: #e8e8e8; --muted: #a0a0a0; --line: #2c2f33;
    --panel: #1c1f23;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: var(--bg); color: var(--fg);
    font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
          Helvetica, Arial, sans-serif;
    padding: 0 16px 48px;
  }}
  .wrap {{ max-width: 980px; margin: 0 auto; }}
  h1 {{ font-size: 1.35rem; margin: 28px 0 4px; font-weight: 650; }}
  .sub {{ color: var(--muted); margin: 0 0 20px; font-size: 0.92rem; }}
  #viewer {{
    width: 100%; height: 60vh; min-height: 380px; position: relative;
    border: 1px solid var(--line); border-radius: 10px; overflow: hidden;
    background: var(--panel);
  }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 14px 0; }}
  button {{
    font: inherit; font-size: 0.86rem; padding: 7px 13px; cursor: pointer;
    background: var(--panel); color: var(--fg);
    border: 1px solid var(--line); border-radius: 7px;
  }}
  button:hover {{ border-color: var(--muted); }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 18px;
           font-size: 0.88rem; }}
  th, td {{ text-align: left; padding: 8px 10px;
            border-bottom: 1px solid var(--line); vertical-align: top; }}
  th {{ color: var(--muted); font-weight: 600; }}
  .swatch {{ display: inline-block; width: 11px; height: 11px;
             border-radius: 3px; margin-right: 7px; vertical-align: -1px; }}
  .caveat {{ margin-top: 26px; padding: 13px 15px; background: var(--panel);
             border: 1px solid var(--line); border-radius: 8px;
             color: var(--muted); font-size: 0.86rem; }}
  code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 0.86em; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>{title}</h1>
  <p class="sub">Drag to rotate, scroll to zoom. Grey cartoon is NS5;
     coloured spheres are the three substitutions that define lineage
     <code>3I_A.2</code>.</p>

  <div id="viewer"></div>

  <div class="controls">
    <button onclick="view('overview')">Overview</button>
    <button onclick="view('interface')">MTase-RdRp interface (H52)</button>
    <button onclick="view('chamber')">Initiation chamber (T404)</button>
    <button onclick="toggleSurface()">Toggle surface</button>
    <button onclick="toggleTheme()">Light / dark</button>
    <button onclick="savePNG()">Download PNG</button>
  </div>
  <p class="sub" style="margin-top:-4px">
    Set the view you want, then Download PNG. Save it as
    <code>results/figures/ns5_structure.png</code> for the README.
  </p>

  <table>
    <thead><tr><th>Substitution</th><th>Domain</th><th>Structural context</th></tr></thead>
    <tbody id="legend"></tbody>
  </table>

  <div class="caveat">
    <strong>What this does and does not show.</strong> Positions and distances
    are measured from the crystal structure. They show where these residues
    sit. They do not show that any of them changes capping, replication or
    fitness &mdash; that would require an enzyme assay, and none was done.
    The initiation chamber exists in the pre-initiation conformation; the
    priming loop retracts during elongation.
  </div>
</div>

<script>
const PDB = {pdb_js};
const HL = {highlight_js};
let viewer, surfaceOn = false, surf = null;

function init() {{
  viewer = $3Dmol.createViewer(document.getElementById("viewer"),
                               {{backgroundColor: "white"}});
  viewer.addModel(PDB, "pdb");

  viewer.setStyle({{}}, {{cartoon: {{color: "#d0d0d0", opacity: 0.92}}}});
  viewer.setStyle({{resi: "1-265"}},   {{cartoon: {{color: "#cfd8dc"}}}});
  viewer.setStyle({{resi: "321-900"}}, {{cartoon: {{color: "#b0bec5"}}}});
  viewer.setStyle({{resi: "{PRIMING_LOOP[0]}-{PRIMING_LOOP[1]}"}},
                  {{cartoon: {{color: "#ffa726"}}}});
  viewer.setStyle({{resi: "{MOTIF_C[0]}-{MOTIF_C[1]}"}},
                  {{cartoon: {{color: "#8e24aa"}}, stick: {{radius: 0.18}}}});
  viewer.setStyle({{resi: {json.dumps(TETRAD)}}},
                  {{cartoon: {{color: "#26a69a"}}, stick: {{radius: 0.18}}}});

  // ligands: SAH, the cap, metals - never the buffer
  viewer.addStyle({{resn: ["SAH", "SAM"]}},
                  {{stick: {{radius: 0.2, colorscheme: "yellowCarbon"}}}});
  viewer.addStyle({{resn: ["M7G", "A", "G", "U", "C"]}},
                  {{stick: {{radius: 0.2, colorscheme: "cyanCarbon"}}}});
  viewer.addStyle({{resn: ["ZN", "MG", "MN"]}}, {{sphere: {{radius: 0.8}}}});

  const tbody = document.getElementById("legend");
  HL.forEach(h => {{
    viewer.addStyle({{resi: h.resi}},
                    {{sphere: {{color: h.color, radius: 1.9}}}});
    viewer.addLabel(h.label, {{
      position: {{resi: h.resi}}, backgroundColor: h.color,
      backgroundOpacity: 0.85, fontColor: "white", fontSize: 11
    }}, {{resi: h.resi}});
    const dom = h.resi <= 265 ? "MTase" : (h.resi <= 320 ? "linker" : "RdRp");
    tbody.insertAdjacentHTML("beforeend",
      `<tr><td><span class="swatch" style="background:${{h.color}}"></span>` +
      `<code>NS5:${{h.label}}</code></td><td>${{dom}}</td><td>${{h.note}}</td></tr>`);
  }});

  viewer.zoomTo();
  viewer.render();
}}

function view(which) {{
  if (which === "overview") viewer.zoomTo();
  if (which === "interface") viewer.zoomTo({{resi: "52"}}, 1200);
  if (which === "chamber")   viewer.zoomTo({{resi: "404"}}, 1200);
  viewer.render();
}}

function toggleSurface() {{
  if (surfaceOn) {{ viewer.removeAllSurfaces(); surfaceOn = false; }}
  else {{
    viewer.addSurface($3Dmol.SurfaceType.VDW,
                      {{opacity: 0.55, color: "#cfd8dc"}});
    surfaceOn = true;
  }}
  viewer.render();
}}

function savePNG() {{
  // 3Dmol renders to a canvas, so the current view can be exported directly.
  // No PyMOL, no extra dependency, and the static image matches the live one.
  const uri = viewer.pngURI();
  const a = document.createElement("a");
  a.href = uri;
  a.download = "ns5_structure.png";
  document.body.appendChild(a);
  a.click();
  a.remove();
}}

function toggleTheme() {{
  const r = document.documentElement;
  const dark = r.getAttribute("data-theme") === "dark";
  r.setAttribute("data-theme", dark ? "light" : "dark");
  viewer.setBackgroundColor(dark ? "white" : "#15171a");
  viewer.render();
}}

init();
</script>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    size_kb = out_path.stat().st_size / 1024
    print(f"  wrote {out_path}  ({size_kb:.0f} KB, self-contained)")


def build_pymol_script(pdb_path, out_path):
    lines = [
        "# Publication figure for the NS5 lineage-defining substitutions.",
        "# Run:  pymol -cq scripts/../results/figures/ns5_structure.pml",
        f"load {pdb_path}, ns5",
        "hide everything",
        "bg_color white",
        "set ray_opaque_background, 1",
        "show cartoon, polymer",
        "color grey80, resi 1-265",
        "color grey60, resi 321-900",
        "color orange, resi 782-809",
        "color purple, resi 662-664",
        "color teal, resi 61+146+180+216",
        "show sticks, resn SAH+SAM+M7G and not hydro",
        "color yellow, resn SAH+SAM",
        "color cyan, resn M7G",
        "show spheres, resn ZN+MG",
        "color slate, resn ZN+MG",
        "show sticks, resn A+G+U+C and polymer.nucleic",
        "color cyan, polymer.nucleic",
        "# buffer is never shown",
        "hide everything, resn ACT+GOL+SO4+EDO+HOH",
    ]
    for resi, meta in HIGHLIGHT.items():
        lines += [
            f"show spheres, resi {resi} and name CA",
            f"set sphere_scale, 0.9, resi {resi} and name CA",
            f"color {meta['color'].replace('#', '0x')}, resi {resi}",
        ]
    lines += [
        "set cartoon_transparency, 0.15",
        "orient polymer",
        "set ray_trace_mode, 0",
        "ray 2000, 1500",
        "png results/figures/ns5_structure.png, dpi=300",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out_path}")


def build_conservation_panel(out_path):
    tsv = Path("results/structure/position_conservation.tsv")
    if not tsv.exists():
        print(f"  skipped conservation panel: {tsv} not found "
              f"(run scripts/19 first)")
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  skipped conservation panel: matplotlib not installed")
        return

    d = pd.read_csv(tsv, sep="\t")
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    labels, pcts, colors = [], [], []
    for _, row in d.iterrows():
        pos = int(str(row["position"]).split(":")[1])
        meta = HIGHLIGHT.get(pos, {"sub": row["position"], "color": "#888888"})
        labels.append(f"NS5:{meta['sub']}\n(n={int(row['n_covered']):,} covered)")
        pcts.append(row["pct_reference_state"])
        colors.append(meta["color"])
    bars = ax.bar(labels, pcts, color=colors, width=0.55)
    for b, p in zip(bars, pcts):
        ax.text(b.get_x() + b.get_width() / 2, p + 0.12, f"{p:.2f}%",
                ha="center", fontsize=9)
    ax.set_ylim(95, 100.6)
    ax.set_ylabel("% of covered DENV-3 sequences\ncarrying the reference residue")
    ax.set_title("All three positions are near-invariant across DENV-3;\n"
                 "lineage 3I_A.2 carries the minority state at each",
                 fontsize=10.5, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="data/structures/5DTO.pdb")
    ap.add_argument("--html", default="docs/ns5_3d.html")
    args = ap.parse_args()

    pdb = Path(args.pdb)
    if not pdb.exists():
        raise SystemExit(f"Structure not found: {pdb}")

    title = (f"NS5 lineage-defining substitutions of DENV-3 3I_A.2 "
             f"({pdb.stem})")
    print("Building figures ...")
    build_3dmol_html(pdb, Path(args.html), title)
    build_pymol_script(pdb, Path("results/figures/ns5_structure.pml"))
    build_conservation_panel(Path("results/figures/ns5_conservation.png"))

    print("\nNext:")
    print("  pymol -cq results/figures/ns5_structure.pml   # static PNG")
    print("  Enable GitHub Pages on /docs, then the viewer is at")
    print("  https://<user>.github.io/<repo>/ns5_3d.html")


if __name__ == "__main__":
    main()
