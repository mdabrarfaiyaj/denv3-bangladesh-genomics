# Publication figure for the NS5 lineage-defining substitutions.
# Run:  pymol -cq scripts/../results/figures/ns5_structure.pml
load data/structures/5DTO.pdb, ns5
hide everything
bg_color white
set ray_opaque_background, 1
show cartoon, polymer
color grey80, resi 1-265
color grey60, resi 321-900
color orange, resi 782-809
color purple, resi 662-664
color teal, resi 61+146+180+216
show sticks, resn SAH+SAM+M7G and not hydro
color yellow, resn SAH+SAM
color cyan, resn M7G
show spheres, resn ZN+MG
color slate, resn ZN+MG
show sticks, resn A+G+U+C and polymer.nucleic
color cyan, polymer.nucleic
# buffer is never shown
hide everything, resn ACT+GOL+SO4+EDO+HOH
show spheres, resi 52 and name CA
set sphere_scale, 0.9, resi 52 and name CA
color 0xe6194b, resi 52
show spheres, resi 139 and name CA
set sphere_scale, 0.9, resi 139 and name CA
color 0x3cb44b, resi 139
show spheres, resi 404 and name CA
set sphere_scale, 0.9, resi 404 and name CA
color 0x4363d8, resi 404
set cartoon_transparency, 0.15
orient polymer
set ray_trace_mode, 0
ray 2000, 1500
png results/figures/ns5_structure.png, dpi=300
