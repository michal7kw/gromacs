# Bugs fixed

## Collective variables module (Colvars) update

The ([Colvars](https://colvars.github.io)) library for enhanced sampling
simulations included in GROMACS has been
updated to version 2025-10-13.

This update brings many fixes, including some errors in gradient
calculations, mostly affecting rarely used cases (see [this page for
details](https://github.com/Colvars/colvars/wiki/Gradient-evaluation-bugs-fixed-in-September-2025)).
A complete list of changes can be found
[here](https://gitlab.com/gromacs/gromacs/-/merge_requests/5397).

## Tpr-file versioning fixed

Older installed GROMACS versions are intended
to be able to do at least basic reading of .tpr files written by newer
versions of GROMACS, unless made impossible by
changes to the .tpr format. This happened in
GROMACS 2025, but was not accompanied an
appropriate increment of the version numbers, resulting in incorrect and
misleading error messages. This cannot be fixed retroactively for
installed GROMACS versions and extant .tpr
files, but at least the problem conditions will no longer continue to
arise for newly created .tpr files.

`5334`

## Fixed missing non-bonded interactions with direct halo communication

When using the experimental direct halo communication feature combined
with 8-wide SIMD non-bonded kernels and OpenMP threading, non-bonded
interactions could be missing.

`5509`

## Allow atoms involved intermolecular-exclusion to be perturbed

`5527`

## Added check for constructing atoms of virtual sites

Constructing atoms for virtual sites can themselves be virtual sites,
but only when those constructing atoms are virtual sites are higher up
in the function type list (i.e. simpler constructions). This was
documented in the manual. Now`grompp` and `mdrun` will throw an error
when these restrictions are violated.

`5535`
