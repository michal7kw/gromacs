# General properties

`gmx energy `, `gmx traj `  
To analyze some or all *energies* and other properties, such as *total
pressure*, *pressure tensor*, *density*, *box-volume* and *box-sizes*,
use the program `gmx energy `. A choice can be made from a
list a set of energies, like potential, kinetic or total energy, or
individual contributions, like Lennard-Jones or dihedral energies.

The *center-of-mass velocity*, defined as

$$
{\bf v}_{com} = {1 \over M} \sum_{i=1}^N m_i {\bf v}_i
$$

with $M = \sum_{i=1}^N m_i$ the total mass of the system, can be
monitored in time by the program `gmx traj ` `-com -ov`. It is
however recommended to remove the center-of-mass velocity every step
(see chapter `algorithms`)!
