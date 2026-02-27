# Improvements to GROMACS tools

## Fixed bug in gmx order -calcdist

The reference position for the distance calculation was calculated
wrongly.

## Improved grompp usability by rejecting more invalid .mdp lines

Lines like

> ref-t 298 = 0.1 =

are now all rejected with a descriptive message, which will help prevent
some kinds of errors in constructing .mdp inputs. Note that an .mdp
parameter name with a missing value is still accepted, and leads to the
default behavior for that parameter.

## Added convert-trj

A new tool `convert-trj ` has been added to allow users
to interchange trajectory formats without having to use legacy
`gmx trjconv`. Supported actions include the generation of slimmed down
output trajectories, as well as the replacement of particle information
in individual frames with data from a structure file. The new tool
allows the usage of command line selections, meaning it is no longer
necessary to write `index <ndx>` files to select certain atoms. It is
part of the drive to split up the `trjconv ` tool into
smaller parts.

## Added extract-cluster

Added a dedicated tool to extract trajectory frames corresponding to
different clusters obtained from `gmx cluster`. The new
`extract-cluster ` tool generates new trajectories
that contain only those frames that correspond to the correct cluster.
The corresponding option **-sub** in `gmx trjconv` has been removed.

## Changed behaviour of genion

Functionality of genion was altered to prevent swapping ions for solvent
closer than -rmin from any other non-solvent atom. This improvement
prevents situations where an ion could be placed at the core of a
protein, which would potentially render the folded protein less stable
or may require long equilibration times.
