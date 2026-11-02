load /home/abrar/projects/denv3-bangladesh-genomics/data/structures/1UZG.pdb
hide everything
show cartoon
color grey80
bg_color white
select subs, chain A and resi 139+219+270+231+301+303
show spheres, subs
color red, subs
show sticks, hetatm and not resn HOH
color yellow, hetatm and not resn HOH
orient
