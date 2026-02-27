# Miscellaneous

## grompp no longer modifies nstcomm

grompp will no longer set nstcomm, the interval for center of mass
motion removal, equal to nstcalcenergy when nstcomm \< nstcalcenergy. A
note is still printed in that case.

## Bonded atom types names can now start with a digit

Bonded atom types names in topologies were not allowed to start with a
number. Now all names are supported that contain at least one non-digit
character.

`4120`

## grompp now warns when exclusion forces might be missing

When using PME, exclusions between non-perturbed,atom pairs should be
within the cut-off distance, otherwise mdrun might not compute grid
correction forces and energies. grompp now computes these distance for
the starting structure and warns when they are beyond 90% of the cut-off
distance and generates an error when they are beyond the cut-off
distance.

`4051`

## The AWH cover diameter for angles now has units degrees

Using old tpr files that apply AWH to angles or dihedrals and have a
non-zero cover diameter results in an error with the suggestion to
regenerate the tpr file.

`4367`

## Core spin-up code is removed

Formerly, on non-x86 and non-PowerPC platforms, mdrun ran some
multi-threaded code to try to wake up any cores that the OS might have
powered down. This caused problems on some Arm platforms, and does not
seem to suit a significant number of platforms for use of
GROMACS. So now it is removed.

If required, please manually spin-up the cores with, e.g.,
`stress --cpu $(nproc --all)`.

`4074`

### Add documentation for linear angle potential

Added documentation and reference for the linear angle potential. Also
added please_cite entry, but there is no call to reference it yet.

`4286`

## gmxapi.mdrun guarantees trajectory output

gmxapi simulations now always run with full-precision trajectory output
(`-o`) in order to guarantee the availability of a usable output
trajectory through the `mdrun.output.trajectory` result.

`4285`

## gmxapi.mdrun accepts arbitrary runtime arguments

Arbitrary mdrun arguments can be passed through gmxapi with the new
*runtime_args* key word argument, accepting a dictionary of flags and
values.

`4284`

## Improved MPI awareness and task uniqueness for gmxapi Python runner

Previously, only the Python components in `gmxapi.simulation` reacted to
the presence of an MPI context. This could result in duplicate work or
even invalid file access.

`gmxapi.commandline_operation` now executes tasks in unique working
directories.

For all gmxapi operations, tasks are only launched from one process (per
ensemble member). If [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
is available, the MPI environment is inspected. If multiple ranks are
discovered, the `ResourceManager` instances on the various ranks
coordinate to make sure that `update` is only called for each member of
each task once. Results are broadcast to all ranks from the
`ResourceManager` where the work occurred.

These changes merely constitute a bug-fix. Additional development is
needed for more optimal use of resources and to reduce unnecessary data
transfers.

`3138`

## Further discouraged use of Berendsen coupling algorithms

Those algorithms have been proven to cause incorrect sampling of their
respective distributions and are mainly provided as a means to provide
backwards compatibility for older simulations. This is why their use has
been further discouraged by changing the current notes about their use
to actual warnings at grompp time.
