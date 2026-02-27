# Protein-related items

`gmx dssp `, `gmx rama `, `gmx wheel `  
To analyze structural changes of a protein, you can calculate the radius
of gyration or the minimum residue distances over time (see sec. `rg`),
or calculate the RMSD (sec. `rmsd`).

To analyze the secondary structure of a protein (not only for static
structures, but also for trajectories), you can use the program
`gmx dssp `, which is a native implementation of DSSP
algorithm `176 <refKabsch83>`, but also is based on the ([DSSP V.4
algorithm](https://github.com/PDB-REDO/dssp)) with some additional
features. For example, you can take into account native hydrogens from
the structure (`-hmode gromacs`, set by default), while in the original
algorithm, hydrogen atoms are set as pseudo-atoms with coordinates based
on the coordinates of the MainChain atoms (`-hmode dssp`). Also, it is
possible to conduct a fast search for neighboring residues using
Neighbor Search (`-nb`, default), instead of the slow enumeration of
protein residues among themselves according to the "each with each"
principle implemented in the original algorithm (`-nonb`).

One other important analysis of proteins is the so-called *Ramachandran
plot*. This is the projection of the structure on the two dihedral
angles $\phi$ and $\psi$ of the protein backbone, see
`Fig. %s <fig-phipsi>`: 

![plots/phipsi.*](plots/phipsi.*)

*Definition of the dihedral angles $ϕ$ and $ψ$ of the protein backbone.*

To evaluate this Ramachandran plot you can use the program
`gmx rama `. A typical output is given in
`Fig. %s <fig-rama>`.

![plots/rama.*%C2%A0](plots/rama.*%C2%A0)

*Ramachandran plot of a small protein.*

When studying $\alpha$-helices it is useful to have a *helical wheel*
projection of your peptide, to see whether a peptide is amphipathic.
This can be done using the `gmx wheel ` program. Two examples
are plotted in `Fig. %s <fig-hprwheel>`.

![plots/hpr-wheel.*](plots/hpr-wheel.*)

*Helical wheel projection of the N-terminal helix of HPr.*
