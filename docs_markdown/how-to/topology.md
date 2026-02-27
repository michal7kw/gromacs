# Adding a Residue to a Force Field

## Adding a new residue

If you have the need to introduce a new residue into an existing force
field so that you can use `pdb2gmx `, or modify an existing
one, there are several files you will need to modify. You must consult
the `/reference-manual/index/` for description of the required format.
Follow these steps:

1.  Add the residue to the `rtp` file for your chosen force field. You
    might be able to copy an existing residue, rename it and modify it
    suitably, or you may need to use an external topology generation
    tool and adapt the results to the `rtp` format.
2.  If you need hydrogens to be able to be added to your residue, create
    an entry in the relevant `hdb` file.
3.  If you are introducing new atom types, add them to the
    `atomtypes.atp` and `ffnonbonded.itp` files.
4.  If you require any new bonded types, add them to `ffbonded.itp`.
5.  Add your residue to `residuetypes.dat` with the appropriate
    specification (Protein, DNA, Ion, etc).
6.  If the residue involves special connectivity to other residues,
    update `specbond.dat`.

Note that if all you are doing is simulating some weird ligand in water,
or some weird ligand with a normal protein, then the above is more work
than generating a standalone `itp` file containing a `[moleculetype]`
(for example, by modifying the `top` produced by some parameterization
server), and inserting an `#include` of that `itp` file into a `top`
generated for the system without that weird ligand.

## Modifying a force field

Modifying a force field is best done by making a full copy of the
installed forcefield directory and `residuetypes.dat` into your local
working directory:

    cp -r $GMXLIB/residuetypes.dat $GMXLIB/amber99sb.ff .

Then, modify those local copies as above. `pdb2gmx ` will
then find both the original and modified version and you can choose the
modified version interactively from the list, or if you use the
`pdb2gmx ` `-ff` option the local version will override the
system version.

# Water solvation

When using `solvate ` to generate a box of solvent, you
need to supply a pre-equilibrated box of a suitable solvent for
`solvate ` to stack around your solute(s), and then to
truncate to give the simulation volume you desire. When using any
3-point model (e.g. `SPC`, `SPC/E` or `TIP3P`) you should specify
`-cs spc216.gro` which will take this file from `the gromacs/share/top`
directory. Other water models (e.g. `TIP4P` and `TIP5P`) are available
as well. Check the contents of the `/share/top` subdirectory of your
GROMACS installation. After solvation, you
should then be sure to equilibrate for at least 5-10ps at the desired
temperature. You will need to select the right water model in your `top`
file, either with the `-water` flag to `pdb2gmx `, or by
editing your `top` file appropriately by hand.

For information about how to use solvents other than pure water, please
see `Non-Water Solvation ` or
`Mixed Solvents `.

# Non water solvent

It is possible to use solvents other than water in
GROMACS. The only requirements are that you
have a pre-equilibrated box of whatever solvent you need, and suitable
parameters for this species in a simulation. One can then pass the
solvent box to the -cs switch of `solvate ` to accomplish
solvation.

A series of about 150 different equilibrated liquids validated for use
with GROMACS, and for the OPLS/AA and GAFF
force fields, can be found at
[virtualchemistry](https://virtualchemistry.org/).

## Making a non-aqueous solvent box

Choose a box density and box size. The size does not have to be that of
your eventual simulation box - a 1nm cube is probably fine. Generate a
single molecule of the solvent. Work out how much volume a single
molecule would have in the box of your chosen density and size. Use
`editconf ` to place a box of that size around your single
molecule. Then use `editconf ` to move the molecule a
little bit off center. Then use `genconf ` `-rot` to
replicate that box into a large one of the right size and density. Then
equilibrate thoroughly to remove the residual ordering of the molecules,
using NVT and periodic boundary conditions. Now you have a box you can
pass to `solvate ` `-cs`, which will replicate it to fit
the size of the actual simulation box.

# Mixed solvent

A common question that new users have is how to create a system with
mixed solvent (urea or DMSO at a given concentration in water, for
example). The simplest procedure for accomplishing this task is as
follows:

- Determine the number of co-solvent molecules necessary, given the box
  dimensions of your system.
- Generate a coordinate file of a single molecule of your co-solvent
  (i.e., `urea.gro`).
- Use the `-ci -nmol` options of `gmx insert-molecules` to add the
  required number of co-solvent molecules to the box.
- Fill the remainder of the box with water (or whatever your other
  solvent is) using `gmx solvate` or `gmx insert-molecules`.
- Edit your `topology ` to `#include` the appropriate `itp` files,
  as well as make changes to the `[ molecules ]` directive to account
  for all the species in your system.

# Making Disulfide Bonds

The easiest way to do this is by using the mechanism implemented with
the `specbond.dat` file and `pdb2gmx `. You may find
`pdb2gmx ` `-ss yes` is useful. The sulfur atoms will need
to be in the same unit that `pdb2gmx ` is converting to a
`moleculetype`, so invoking `pdb2gmx ` `-chainsep`
correctly may be required. See `pdb2gmx ` `-h`. This
requires that the two sulfur atoms be within a distance + tolerance
(usually 10%) in order to be recognised as a disulfide. If your sulfur
atoms are not this close, then either you can

- edit the contents of `specbond.dat` to allow the bond formation and do
  energy minimization very carefully to allow the bond to relax to a
  sensible length, or
- run a preliminary EM or MD with a distance restraint (and no disulfide
  bond) between these sulfur atoms with a large force constant so that
  they approach within the existing `specbond.dat` range to provide a
  suitable coordinate file for a second invocation of
  `pdb2gmx `.

Otherwise, editing your `top` file by hand is the only option.
