# `pdb2gmx ` input files

The GROMACS program `pdb2gmx `
generates a topology for the input coordinate file. Several formats are
supported for that coordinate file, but `pdb` is the most commonly-used
format (hence the name `pdb2gmx `). `pdb2gmx `
searches for force fields in sub-directories of the
GROMACS `share/top` directory and your working
directory. Force fields are recognized from the file `forcefield.itp` in
a directory with the extension `.ff`. The file `forcefield.doc` may be
present, and if so, its first line will be used by
`pdb2gmx ` to present a short description to the user to
help in choosing a force field. Otherwise, the user can choose a force
field with the `-ff xxx` command-line argument to
`pdb2gmx `, which indicates that a force field in a
`xxx.ff` directory is desired. `pdb2gmx ` will search first
in the working directory, then in the GROMACS
`share/top` directory, and use the first matching `xxx.ff` directory
found.

Two general files are read by `pdb2gmx `: an atom type file
(extension `atp`, see `atomtype`) from the force-field directory, and a
file called `residuetypes.dat` from either the working directory, or the
GROMACS `share/top` directory.
`residuetypes.dat` determines which residue names are considered
protein, DNA, RNA, water, and ions.

`pdb2gmx ` can read one or multiple databases with
topological information for different types of molecules. A set of files
belonging to one database should have the same basename, preferably
telling something about the type of molecules (*e.g.* aminoacids, rna,
dna). The possible files are:

- `<basename>.rtp`
- `<basename>.r2b (optional)`
- `<basename>.arn (optional)`
- `<basename>.hdb (optional)`
- `<basename>.n.tdb (optional)`
- `<basename>.c.tdb (optional)`

Only the `rtp` file, which contains the topologies of the building
blocks, is mandatory. Information from other files will only be used for
building blocks that come from an `rtp` file with the same base name.
The user can add building blocks to a force field by having additional
files with the same base name in their working directory. By default,
only extra building blocks can be defined, but calling
`pdb2gmx ` with the `-rtpo` option will allow building
blocks in a local file to replace the default ones in the force field.

## Residue database

The files holding the residue databases have the extension `rtp`.
Originally this file contained building blocks (amino acids) for
proteins, and is the GROMACS interpretation of
the `rt37c4.dat` file of GROMOS. So the residue database file contains
information (bonds, charges, charge groups, and improper dihedrals) for
a frequently-used building block. It is better *not* to change this file
because it is standard input for `pdb2gmx `, but if changes
are needed make them in the `top` file (see `topfile`), or in a `rtp`
file in the working directory as explained in sec. `pdb2gmxfiles`.
Defining topologies of new small molecules is probably easier by writing
an include topology file `itp` directly. This will be discussed in
section `molitp`. When adding a new protein residue to the database, do
not forget to add the residue name to the residuetypes.dat file, so that
`grompp `, `make_ndx ` and analysis tools can
recognize the residue as a protein residue (see `defaultgroups`).

The `rtp` files are only used by `pdb2gmx `. As mentioned
before, the only extra information this program needs from the `rtp`
database is bonds, charges of atoms, charge groups, and improper
dihedrals, because the rest is read from the coordinate input file. Some
proteins contain residues that are not standard, but are listed in the
coordinate file. You have to construct a building block for this
“strange” residue, otherwise you will not obtain a `top` file. This also
holds for molecules in the coordinate file such as ligands, polyatomic
ions, crystallization co-solvents, etc. The residue database is
constructed in the following way:

    [ bondedtypes ]  ; mandatory
    ; bonds  angles  dihedrals  impropers
         1       1          1          2  ; mandatory

    [ GLY ]  ; mandatory

     [ atoms ]  ; mandatory
    ; name  type  charge  chargegroup
         N     N  -0.280     0
         H     H   0.280     0
        CA   CH2   0.000     1
         C     C   0.380     2
         O     O  -0.380     2

     [ bonds ]  ; optional
    ;atom1 atom2      b0      kb
         N     H
         N    CA
        CA     C
         C     O
        -C     N

     [ exclusions ]  ; optional
    ;atom1 atom2

     [ angles ]  ; optional
    ;atom1 atom2 atom3    th0    cth

     [ dihedrals ]  ; optional
    ;atom1 atom2 atom3 atom4   phi0     cp   mult

     [ impropers ]  ; optional
    ;atom1 atom2 atom3 atom4     q0     cq
         N    -C    CA     H
        -C   -CA     N    -O

    [ ZN ]

     [ atoms ]
        ZN    ZN   2.000     0

The file is free format; the only restriction is that there can be at
most one entry on a line. The first field in the file is the
`[ bondedtypes ]` field, which is followed by four numbers, indicating
the interaction type for bonds, angles, dihedrals, and improper
dihedrals. The file contains residue entries, which consist of atoms and
(optionally) bonds, angles, dihedrals, and impropers. The charge group
codes denote the charge group numbers. Atoms in the same charge group
should always be ordered consecutively. When using the hydrogen database
with `pdb2gmx ` for adding missing hydrogens (see `hdb`),
the atom names defined in the `rtp` entry should correspond exactly to
the naming convention used in the hydrogen database. The atom names in
the bonded interaction can be preceded by a minus or a plus, indicating
that the atom is in the preceding or following residue respectively.
Explicit parameters added to bonds, angles, dihedrals, and impropers
override the standard parameters in the `itp` files. This should only be
used in special cases. Instead of parameters, a string can be added for
each bonded interaction. This is used in GROMOS-96 `rtp` files. These
strings are copied to the topology file and can be replaced by
force-field parameters by the C-preprocessor in `grompp `
using `#define` statements.

`pdb2gmx ` automatically generates all angles. This means
that for most force fields the `[ angles ]` field is only useful for
overriding `itp` parameters. For the GROMOS-96 force field the
interaction number of all angles needs to be specified.

`pdb2gmx ` automatically generates one proper dihedral for
every rotatable bond, preferably on heavy atoms. When the
`[ dihedrals ]` field is used, no other dihedrals will be generated for
the bonds corresponding to the specified dihedrals. It is possible to
put more than one dihedral function on a rotatable bond. In the case of
CHARMM27 FF `pdb2gmx ` can add correction maps to the
dihedrals using the default `-cmap` option. Please refer to `charmmff`
for more information.

`pdb2gmx ` sets the number of exclusions to 3, which means
that interactions between atoms connected by at most 3 bonds are
excluded. Pair interactions are generated for all pairs of atoms that
are separated by 3 bonds (except pairs of hydrogens). When more
interactions need to be excluded, or some pair interactions should not
be generated, an `[ exclusions ]` field can be added, followed by pairs
of atom names on separate lines. All non-bonded and pair interactions
between these atoms will be excluded.

## Residue to building block database

Each force field has its own naming convention for residues. Most
residues have consistent naming, but some, especially those with
different protonation states, can have many different names. The `r2b`
files are used to convert standard residue names to the force-field
build block names. If no `r2b` is present in the force-field directory
or a residue is not listed, the building block name is assumed to be
identical to the residue name. The `r2b` can contain 2 or 5 columns. The
2-column format has the residue name in the first column and the
building block name in the second. The 5-column format has 3 additional
columns with the building block for the residue occurring in the
N-terminus, C-terminus and both termini at the same time (single residue
molecule). This is useful for, for instance, the AMBER force fields. If
one or more of the terminal versions are not present, a dash should be
entered in the corresponding column.

There is a GROMACS naming convention for
residues which is only apparent (except for the `pdb2gmx `
code) through the `r2b` file and `specbond.dat` files. This convention
is only of importance when you are adding residue types to an `rtp`
file. The convention is listed in `Table %s <tab-r2b>`. For special
bonds with, for instance, a heme group, the
GROMACS naming convention is introduced
through `specbond.dat` (see `specbond`), which can subsequently be
translated by the `r2b` file, if required.


Internal GROMACS residue
naming convention.


GROMACS ID
Residue


ARG
protonated arginine


ARGN
neutral arginine


ASP
negatively charged aspartic acid


ASPH
neutral aspartic acid


CYS
neutral cysteine


CYS2
cysteine with sulfur bound to another cysteine or a
heme


GLU

negatively charged glutamic acid


GLUH

neutral glutamic acid


HISD
neutral histidine with N$_{*δ*}$ protonated


HISE
neutral histidine with N$_{*ϵ*}$ protonated


HISH
positive histidine with both N$_{*δ*}$ and N$_{*ϵ*}$ protonated


HIS1
histidine bound to a heme


LYSN
neutral lysine


LYS
protonated lysine


HEME
heme


## Atom renaming database

Force fields often use atom names that do not follow IUPAC or PDB
convention. The `arn` database is used to translate the atom names in
the coordinate file to the force-field names. Atoms that are not listed
keep their names. The file has three columns: the building block name,
the old atom name, and the new atom name, respectively. The residue name
supports question-mark wildcards that match a single character.

An additional general atom renaming file called `xlateat.dat` is present
in the `share/top` directory, which translates common non-standard atom
names in the coordinate file to IUPAC/PDB convention. Thus, when writing
force-field files, you can assume standard atom names and no further
atom name translation is required, except for translating from standard
atom names to the force-field ones.

## Hydrogen database

The hydrogen database is stored in `hdb` files. It contains information
for the `pdb2gmx ` program on how to connect hydrogen atoms
to existing atoms. In versions of the database before
GROMACS 3.3, hydrogen atoms were named after
the atom they are connected to: the first letter of the atom name was
replaced by an ‘H.’ In the versions from 3.3 onwards, the H atom has to
be listed explicitly, because the old behavior was protein-specific and
hence could not be generalized to other molecules. If more than one
hydrogen atom is connected to the same atom, a number will be added to
the end of the hydrogen atom name. For example, adding two hydrogen
atoms to `ND2` (in asparagine), the hydrogen atoms will be named `HD21`
and `HD22`. This is important since atom naming in the `rtp` file must
be the same. The format of the hydrogen database is as follows:

    ; res   # additions
            # H add type    H       i       j       k
    ALA     1
            1       1       H       N       -C      CA
    ARG     4
            1       2       H       N       CA      C
            1       1       HE      NE      CD      CZ
            2       3       HH1     NH1     CZ      NE
            2       3       HH2     NH2     CZ      NE

On the first line we see the residue name (ALA or ARG) and the number of
kinds of hydrogen atoms that may be added to this residue by the
hydrogen database. After that follows one line for each addition, on
which we see:

- The number of H atoms added
- The method for adding H atoms, which can be any of:
  1.  *one planar hydrogen*, *e.g.* *rings or peptide bond*  
      One hydrogen atom (n) is generated, lying in the plane of atoms
      (i,j,k) on the plane bisecting angle (j-i-k) at a distance of 0.1
      nm from atom i, such that the angles (n-i-j) and (n-i-k) are $>$
      90$^{\rm o}$.

  2.  *one single hydrogen*, *e.g.* *hydroxyl*  
      One hydrogen atom (n) is generated at a distance of 0.1 nm from
      atom i, such that angle (n-i-j)=109.5 degrees and dihedral
      (n-i-j-k)=trans.

  3.  *two planar hydrogens*, *e.g.* *ethylene -C=CH*$_2$, *or amide
      -C(=O)NH*$_2$  
      Two hydrogens (n1,n2) are generated at a distance of 0.1 nm from
      atom i, such that angle (n1-i-j)=(n2-i-j)=120 degrees and dihedral
      (n1-i-j-k)=cis and (n2-i-j-k)=trans, such that names are according
      to IUPAC standards `129 <refiupac70>`.

  4.  *two or three tetrahedral hydrogens*, *e.g.* *-CH*$_3$  
      Three (n1,n2,n3) or two (n1,n2) hydrogens are generated at a
      distance of 0.1 nm from atom i, such that angle
      (n1-i-j)=(n2-i-j)=(n3-i-j)=109.47$^{\rm o}$, dihedral
      (n1-i-j-k)=trans, (n2-i-j-k)=trans+120 and
      (n3-i-j-k)=trans+240$^{\rm o}$.

  5.  *one tetrahedral hydrogen*, *e.g.* *C*$_3$\*CH\*  
      One hydrogen atom (n$^\prime$) is generated at a distance of 0.1
      nm from atom i in tetrahedral conformation such that angle
      (n$^\prime$-i-j)=(n$^\prime$-i-k)=(n$^\prime$-i-l)=109.47$^{\rm o}$.

  6.  *two tetrahedral hydrogens*, *e.g.* *C-CH*$_2$\*-C\*  
      Two hydrogen atoms (n1,n2) are generated at a distance of 0.1 nm
      from atom i in tetrahedral conformation on the plane bisecting
      angle j-i-k with angle
      (n1-i-n2)=(n1-i-j)=(n1-i-k)=109.47$^{\rm o}$.

  7.  *two water hydrogens*  
      Two hydrogens are generated around atom i according to
      SPC `80 <refBerendsen81>` water geometry. The symmetry axis will
      alternate between three coordinate axes in both directions.

  8.  *three water “hydrogens”*  
      Two hydrogens are generated around atom i according to
      SPC `80 <refBerendsen81>` water geometry. The symmetry axis will
      alternate between three coordinate axes in both directions. In
      addition, an extra particle is generated on the position of the
      oxygen with the first letter of the name replaced by ‘M’. This is
      for use with four-atom water models such as
      TIP4P `128 <refJorgensen83>`.

  9.  *four water “hydrogens”*  
      Same as above, except that two additional particles are generated
      on the position of the oxygen, with names ‘LP1’ and ‘LP2.’ This is
      for use with five-atom water models such as
      TIP5P `130 <refMahoney2000a>`.
- The name of the new H atom (or its prefix, *e.g.* `HD2` for the
  asparagine example given earlier).
- Three or four control atoms (i,j,k,l), where the first always is the
  atom to which the H atoms are connected. The other two or three depend
  on the code selected. For water, there is only one control atom.

Some more exotic cases can be approximately constructed from the above
tools, and with suitable use of energy minimization are good enough for
beginning MD simulations. For example secondary amine hydrogen, nitrenyl
hydrogen ($\mathrm{C}=\mathrm{NH}$) and even ethynyl hydrogen could be
approximately constructed using method 2 above for hydroxyl hydrogen.

## Termini database

The termini databases are stored in `aminoacids.n.tdb` and
`aminoacids.c.tdb` for the N- and C-termini respectively. They contain
information for the `pdb2gmx ` program on how to connect
new atoms to existing ones, which atoms should be removed or changed,
and which bonded interactions should be added. Their format is as
follows (from `gromos43a1.ff/aminoacids.c.tdb`):

    [ None ]

    [ COO- ]
    [ replace ]
    C   C   C   12.011  0.27
    O   O1  OM  15.9994 -0.635
    OXT O2  OM  15.9994 -0.635
    [ add ]
    2   8   O   C   CA  N
        OM  15.9994 -0.635
    [ bonds ]
    C   O1  gb_5
    C   O2  gb_5
    [ angles ]
    O1  C   O2  ga_37
    CA  C   O1  ga_21
    CA  C   O2  ga_21
    [ dihedrals ]
    N   CA  C   O2  gd_20
    [ impropers ]
    C   CA  O2  O1  gi_1

The file is organized in blocks, each with a header specifying the name
of the block. These blocks correspond to different types of termini that
can be added to a molecule. In this example `[ COO- ]` is the first
block, corresponding to changing the terminal carbon atom into a
deprotonated carboxyl group. `[ None ]` is the second terminus type,
corresponding to a terminus that leaves the molecule as it is. Block
names cannot be any of the following: `replace`, `add`, `delete`,
`bonds`, `angles`, `dihedrals`, `impropers`. Doing so would interfere
with the parameters of the block, and would probably also be very
confusing to human readers.

For each block the following options are present:

- `[ replace ]`  
  Replace an existing atom by one with a different atom type, atom name,
  charge, and/or mass. This entry can be used to replace an atom that is
  present both in the input coordinates and in the `rtp` database, but
  also to only rename an atom in the input coordinates such that it
  matches the name in the force field. In the latter case, there should
  also be a corresponding `[ add ]` section present that gives
  instructions to add the same atom, such that the position in the
  sequence and the bonding is known. Such an atom can be present in the
  input coordinates and kept, or not present and constructed by
  `pdb2gmx `. For each atom to be replaced on line should
  be entered with the following fields:

  - name of the atom to be replaced
  - new atom name (optional)
  - new atom type
  - new mass
  - new charge

- `[ add ]`  
  Add new atoms. For each (group of) added atom(s), a two-line entry is
  necessary. The first line contains the same fields as an entry in the
  hydrogen database (name of the new atom, number of atoms, type of
  addition, control atoms, see `hdb`), but the possible types of
  addition are extended by two more, specifically for C-terminal
  additions:

  1.  *two carboxyl oxygens, -COO*$^-$  
      Two oxygens (n1,n2) are generated according to rule 3, at a
      distance of 0.136 nm from atom i and an angle
      (n1-i-j)=(n2-i-j)=117 degrees

  2.  *carboxyl oxygens and hydrogen, -COOH*  
      Two oxygens (n1,n2) are generated according to rule 3, at
      distances of 0.123 nm and 0.125 nm from atom i for n1 and n2,
      respectively, and angles (n1-i-j)=121 and (n2-i-j)=115 degrees.
      One hydrogen (n$^\prime$) is generated around n2 according to
      rule 2, where n-i-j and n-i-j-k should be read as
      n$^\prime$-n2-i and n$^\prime$-n2-i-j, respectively.

  After this line, another line follows that specifies the details of
  the added atom(s), in the same way as for replacing atoms, *i.e.*:

  - atom type
  - mass
  - charge
  - charge group (optional)

  Like in the hydrogen database (see `rtp`), when more than one atom is
  connected to an existing one, a number will be appended to the end of
  the atom name. **Note** that, like in the hydrogen database, the atom
  name is now on the same line as the control atoms, whereas it was at
  the beginning of the second line prior to
  GROMACS version 3.3. When the charge group
  field is left out, the added atom will have the same charge group
  number as the atom that it is bonded to.

- `[ delete ]`  
  Delete existing atoms. One atom name per line.

- `[ bonds ]`, `[ angles ]`, `[ dihedrals ]` and `[ impropers ]`  
  Add additional bonded parameters. The format is identical to that used
  in the `rtp` file, see `rtp`.

## Virtual site database

Since we cannot rely on the positions of hydrogens in input files, we
need a special input file to decide the geometries and parameters with
which to add virtual site hydrogens. For more complex virtual site
constructs (*e.g.* when entire aromatic side chains are made rigid) we
also need information about the equilibrium bond lengths and angles for
all atoms in the side chain. This information is specified in the `vsd`
file for each force field. Just as for the termini, there is one such
file for each class of residues in the `rtp` file.

The virtual site database is a simple list of information. The first
couple of sections specify which mass centers (typically called
MCH$_3$/MNH$_3$) to use for CH$_3$, NH$_3$, and NH$_2$ groups.
Depending on the equilibrium bond lengths and angles between the
hydrogens and heavy atoms we need to apply slightly different constraint
distances between these mass centers. **Note** that we do *not* have to
specify the actual parameters (that is automatic), just the type of mass
center to use. To accomplish this, there are three sections names
`[ CH3 ]`, `[ NH3 ]`, and `[ NH2 ]`. For each of these we expect three
columns. The first column is the atom type bound to the 2/3 hydrogens,
the second column is the next heavy atom type which this is bound, and
the third column the type of mass center to use. As a special case, in
the `[ NH2 ]` section it is also possible to specify `planar` in the
second column, which will use a different construction without mass
center. There are currently different opinions in some force fields
whether an NH$_2$ group should be planar or not, but we try hard to
stick to the default equilibrium parameters of the force field.

The second part of the virtual site database contains explicit
equilibrium bond lengths and angles for pairs/triplets of atoms in
aromatic side chains. These entries are currently read by specific
routines in the virtual site generation code, so if you would like to
extend it *e.g.* to nucleic acids you would also need to write new code
there. These sections are named after the short amino acid names
(`[ PHE ]`, `[ TYR ]`, `[ TRP ]`, `[ HID ]`, `[ HIE ]`, `[ HIP ]`), and
simply contain 2 or 3 columns with atom names, followed by a number
specifying the bond length (in nm) or angle (in degrees). **Note** that
these are approximations of the equilibrated geometry for the entire
molecule, which might not be identical to the equilibrium value for a
single bond/angle if the molecule is strained.

## Special bonds

The primary mechanism used by `pdb2gmx ` to generate
inter-residue bonds relies on head-to-tail linking of backbone atoms in
different residues to build a macromolecule. In some cases (*e.g.*
disulfide bonds, a heme group, branched polymers), it is necessary to
create inter-residue bonds that do not lie on the backbone. The file
`specbond.dat` takes care of this function. It is necessary that the
residues belong to the same `[ moleculetype ]`. The `-merge` and
`-chainsep` functions of `pdb2gmx ` can be useful when
managing special inter-residue bonds between different chains.

The first line of `specbond.dat` indicates the number of entries that
are in the file. If you add a new entry, be sure to increment this
number. The remaining lines in the file provide the specifications for
creating bonds. For these bonds, you can also optionally specify a
custom improper dihedral associated with the new bond. The format of the
lines, with optional entries in \[\], is as follows:

`resA atomA nbondsA resB atomB nbondsB length newresA newresB [atomI atomJ atomK atomL]`

The columns indicate:

1.  `resA` The name of residue A that participates in the bond.

2.  `atomA` The name of the atom in residue A that forms the bond.

3.  `nbondsA` The total number of bonds `atomA` can form.

4.  `resB` The name of residue B that participates in the bond.

5.  `atomB` The name of the atom in residue B that forms the bond.

6.  `nbondsB` The total number of bonds `atomB` can form.

7.  `length` The reference length for the bond. If `atomA` and `atomB`
    are not within `length` $\pm$ 10% in the coordinate file supplied
    to `pdb2gmx `, no bond will be formed.

8.  `newresA` The new name of residue A, if necessary. Some force fields
    use *e.g.* CYS2 for a cysteine in a disulfide or heme linkage.

9.  `newresB` The new name of residue B, likewise.

10. `atomI` Custom improper dihedral atom i of i-j-k-l. Has format  
    \[specbond residue\]-\[atom name\] (e.g. B-SG). The letter is either
    A or B corresponding to resA or resB, respectively.

11. `atomJ` Custom improper dihedral atom j of i-j-k-l.

12. `atomK` Custom improper dihedral atom k of i-j-k-l.

13. `atomL` Custom improper dihedral atom l of i-j-k-l.
