# New and improved features

## mdrun now also reports the conserved energy quantity with AWH bias sharing

## Added option for setting the ensemble temperature

Several algorithms, such as pressure coupling and AWH, need the
temperature of the system. When not all atoms are coupled to the (same)
temperature, it is now possible to tell mdrun what the ensemble
temperature is using two new mdp options.

`3854`

## gmxapi.mdrun now publishes the simulation working directory path

`gmxapi.mdrun().output.directory` provides the (Future) path to the
working directory/directories for the simulation(s). This can be useful
in conjunction with `gmxapi.utility.join_path` to express data flow
based on files known to the user to be produced by the simulation but
not represented by other existing attributes of the
`~gmxapi.simulation.mdrun.OutputDataProxy`.

`4548`

## gmxapi.mdrun now captures STDOUT and STDERR

The GROMACS library prints a lot of output
directly to standard out and standard error. Previously, this meant that
simulator output that traditionally goes to the terminal would have to
be caught from outside the Python interpreter. In mpi4py based
ensembles, it could be challenging to catch the output at all, without
manipulating the `mpiexec` command line.

`gmxapi.mdrun` now redirects STDERR and STDOUT during simulation, and
provides paths to the resulting text files on new *stdout* and *stderr*
outputs.

Reference `4541`
