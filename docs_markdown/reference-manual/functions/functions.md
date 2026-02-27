# Interaction function and force fields

To accommodate the potential functions used in some popular force fields
(see `ff`), GROMACS offers a choice of
functions, both for non-bonded interaction and for dihedral
interactions. They are described in the appropriate subsections.

The potential functions can be subdivided into three parts

1.  *Non-bonded*: Lennard-Jones or Buckingham, and Coulomb or modified
    Coulomb. The non-bonded interactions are computed on the basis of a
    neighbor list (a list of non-bonded atoms within a certain radius),
    in which exclusions are already removed.
2.  *Bonded*: covalent bond-stretching, angle-bending, improper
    dihedrals, and proper dihedrals. These are computed on the basis of
    fixed lists.
3.  *Restraints*: position restraints, angle restraints, distance
    restraints, orientation restraints and dihedral restraints, all
    based on fixed lists.
4.  *Applied Forces*: externally applied forces, see chapter `special`.

