# Bugs fixed

## gmx mdrun -append now requires that a checkpoint is found

Previously `gmx mdrun -append` would start from the .tpr configuration
(and thus not append) when the checkpoint file was missing.

## The Verlet buffer now correctly handles perturbed constraints

With free-energy calculations with perturbed constraints, the Verlet
buffer could be underestimated when constraint lengths were perturbed.
As usually only very few constraints are perturbed, the effect is very
small and much smaller than the overestimate of the buffer due to
approximations, so the results of most runs with perturbed constraints
will not have been affected.

`4395`
