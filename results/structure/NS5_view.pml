load /home/abrar/projects/denv3-bangladesh-genomics/data/structures/4V0Q.pdb
hide everything
show cartoon
color grey80
bg_color white
select subs, chain A and resi 52+139+404
show spheres, subs
color red, subs
show sticks, hetatm and not resn HOH
color yellow, hetatm and not resn HOH
orient
