#!/usr/bin/env python3
"""
22_logo.py
----------
Builds the project logo: docs/logo.svg (+ .png).

Hexagonal sticker in the style of other bioinformatics project logos: gradient
fill, bold wordmark on top, organism in the middle, molecular motif at the
bottom.

Iconography chosen for this project specifically:
  * a mature flavivirus particle, drawn as the smooth icosahedral surface with
    its herringbone arrangement of E-protein dimers. Not a phage: dengue is an
    enveloped flavivirus and the surface pattern is its recognisable feature.
  * a phylogeny underneath, because the answer this project gives is a tree.
  * green and red, after the Bangladesh flag, so the logo says where as well
    as what.

The herringbone is generated from ring geometry rather than drawn by hand, so
the dimer count and angles can be tuned in one place.

Run:
  python scripts/22_logo.py
  python scripts/22_logo.py --size 512 --no-png
"""

import argparse
import math
import shutil
import subprocess
from pathlib import Path

# Bangladesh flag green and red, darkened slightly so white type sits on them.
GREEN_DEEP = "#00382C"
GREEN = "#00694F"
RED = "#F42A41"
RED_DEEP = "#C41230"
MINT = "#8FE3C4"
PALE = "#E8FFF6"


def hexagon(cx, cy, r):
    """Flat-top hexagon: vertices at left and right, flat top and bottom."""
    pts = []
    for ang in (180, 240, 300, 0, 60, 120):
        a = math.radians(ang)
        pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)


def virion(cx, cy, r):
    """A flavivirus particle: envelope plus a herringbone of E-protein dimers.

    Each dimer is a short rounded bar. They are placed on concentric rings and
    alternately tilted, which is what gives the real particle its herringbone
    look without trying to draw all 90 dimers accurately.
    """
    out = [f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="url(#virion)" '
           f'stroke="{PALE}" stroke-width="3" opacity="0.97"/>']
    rings = [(0.30, 6), (0.55, 10), (0.78, 14)]
    for frac, n in rings:
        rr = r * frac
        bar = r * 0.26
        for i in range(n):
            a = 2 * math.pi * i / n + (0.22 if frac > 0.5 else 0)
            x = cx + rr * math.cos(a)
            y = cy + rr * math.sin(a)
            tilt = math.degrees(a) + (34 if i % 2 == 0 else -34)
            out.append(
                f'<rect x="{x - bar/2:.1f}" y="{y - r*0.045:.1f}" '
                f'width="{bar:.1f}" height="{r*0.09:.1f}" rx="{r*0.045:.1f}" '
                f'fill="{PALE}" opacity="0.55" '
                f'transform="rotate({tilt:.1f} {x:.1f} {y:.1f})"/>')
    # centre dimer pair, slightly brighter, to anchor the eye
    out.append(f'<circle cx="{cx}" cy="{cy}" r="{r*0.13:.1f}" '
               f'fill="{PALE}" opacity="0.75"/>')
    return "\n".join(out)


def phylogeny(cx, cy, w, h):
    """A small cladogram: one root, splitting right, tips aligned.

    Four tips, one of them highlighted, which is the shape of the finding:
    many lineages, one that persisted.
    """
    sw = max(3.0, w * 0.018)
    left = cx - w / 2
    right = cx + w / 2
    tips = [cy - h * 0.42, cy - h * 0.14, cy + h * 0.14, cy + h * 0.42]
    n1 = (tips[0] + tips[1]) / 2          # node joining tips 0,1
    n2 = (tips[2] + tips[3]) / 2          # node joining tips 2,3
    xr, x1, x2 = left, left + w * 0.34, left + w * 0.60

    def ln(x1_, y1_, x2_, y2_, col=PALE, op=0.85):
        return (f'<line x1="{x1_:.1f}" y1="{y1_:.1f}" x2="{x2_:.1f}" '
                f'y2="{y2_:.1f}" stroke="{col}" stroke-width="{sw:.1f}" '
                f'stroke-linecap="round" opacity="{op}"/>')

    o = [ln(xr, cy, x1, cy),                       # root
         ln(x1, n1, x1, n2),                       # first split
         ln(x1, n1, x2, n1), ln(x1, n2, x2, n2),
         ln(x2, tips[0], x2, tips[1]),
         ln(x2, tips[2], x2, tips[3]),
         ln(x2, tips[0], right, tips[0]),
         ln(x2, tips[1], right, tips[1]),
         ln(x2, tips[2], right, tips[2])]
    # the persisting lineage, in red, running to the present
    o.append(ln(x2, tips[3], right, tips[3], col=RED, op=1.0))
    o.append(f'<circle cx="{right:.1f}" cy="{tips[3]:.1f}" r="{sw*1.6:.1f}" '
             f'fill="{RED}"/>')
    for t in tips[:3]:
        o.append(f'<circle cx="{right:.1f}" cy="{t:.1f}" r="{sw*1.1:.1f}" '
                 f'fill="{PALE}" opacity="0.85"/>')
    return "\n".join(o)


def build(S=1024):
    c = S / 2
    R = S * 0.459
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{S}" height="{S}" '
         f'viewBox="0 0 {S} {S}">',
         '<defs>',
         f'  <linearGradient id="bg" x1="0" y1="0" x2="0.35" y2="1">'
         f'<stop offset="0%" stop-color="{GREEN_DEEP}"/>'
         f'<stop offset="55%" stop-color="{GREEN}"/>'
         f'<stop offset="100%" stop-color="{RED_DEEP}"/></linearGradient>',
         f'  <linearGradient id="word" x1="0" y1="0" x2="0" y2="1">'
         f'<stop offset="0%" stop-color="#FFFFFF"/>'
         f'<stop offset="100%" stop-color="{MINT}"/></linearGradient>',
         f'  <radialGradient id="virion" cx="0.36" cy="0.30" r="0.80">'
         f'<stop offset="0%" stop-color="#FF7A6B"/>'
         f'<stop offset="60%" stop-color="{RED}"/>'
         f'<stop offset="100%" stop-color="{RED_DEEP}"/></radialGradient>',
         '</defs>']

    # sticker: white keyline, then the hexagon
    p.append(f'<polygon points="{hexagon(c, c, R)}" fill="#FFFFFF" '
             f'stroke="#FFFFFF" stroke-width="{S*0.055:.1f}" '
             f'stroke-linejoin="round"/>')
    p.append(f'<polygon points="{hexagon(c, c, R*0.962)}" fill="url(#bg)" '
             f'stroke="url(#bg)" stroke-width="{S*0.03:.1f}" '
             f'stroke-linejoin="round"/>')

    # wordmark
    p.append(f'<text x="{c}" y="{S*0.255:.0f}" text-anchor="middle" '
             f'font-family="Helvetica, Arial, sans-serif" font-weight="800" '
             f'font-size="{S*0.145:.0f}" letter-spacing="{S*0.004:.1f}" '
             f'fill="url(#word)">DENV-3</text>')
    p.append(f'<text x="{c}" y="{S*0.345:.0f}" text-anchor="middle" '
             f'font-family="Helvetica, Arial, sans-serif" font-weight="700" '
             f'font-size="{S*0.072:.0f}" letter-spacing="{S*0.012:.1f}" '
             f'fill="url(#word)" opacity="0.94">BANGLADESH</text>')

    p.append(virion(c, S * 0.545, S * 0.145))
    p.append(phylogeny(c, S * 0.785, S * 0.46, S * 0.17))

    p.append('</svg>')
    return "\n".join(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/logo.svg")
    ap.add_argument("--size", type=int, default=1024)
    ap.add_argument("--no-png", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(args.size))
    print(f"wrote {out}  ({args.size} x {args.size})")

    if args.no_png:
        return
    png = out.with_suffix(".png")
    for cmd in (["rsvg-convert", "-w", str(args.size), str(out), "-o", str(png)],
                ["inkscape", str(out), "--export-filename", str(png),
                 "--export-width", str(args.size)],
                ["convert", "-density", "200", "-background", "none",
                 str(out), str(png)]):
        if shutil.which(cmd[0]):
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                print(f"wrote {png}  (via {cmd[0]})")
                return
            except subprocess.CalledProcessError as e:
                print(f"  {cmd[0]} failed: {e.stderr.decode()[:200]}")
    try:
        import cairosvg
        cairosvg.svg2png(url=str(out), write_to=str(png),
                         output_width=args.size)
        print(f"wrote {png}  (via cairosvg)")
        return
    except ImportError:
        pass
    print("No SVG converter found. GitHub renders SVG in a README directly, "
          "so the PNG is optional.")


if __name__ == "__main__":
    main()
