# Interface-related items

`gmx order `, `gmx density `,
`gmx potential `, `gmx traj `  
When simulating molecules with long carbon tails, it can be interesting
to calculate their average orientation. There are several flavors of
order parameters, most of which are related. The program
`gmx order ` can calculate order parameters using the
equation:

$$
S_{z} = \frac{3}{2}\langle {\cos^2{\theta_z}} \rangle - \frac{1}{2}
$$

where $\theta_z$ is the angle between the $z$-axis of the simulation
box and the molecular axis under consideration. The latter is defined as
the vector from C$_{n-1}$ to C$_{n+1}$. The parameters $S_x$ and
$S_y$ are defined in the same way. The brackets imply averaging over
time and molecules. Order parameters can vary between 1 (full order
along the interface normal) and $-1/2$ (full order perpendicular to
the normal), with a value of zero in the case of isotropic orientation.

The program can do two things for you. It can calculate the order
parameter for each CH$_2$ segment separately, for any of three axes,
or it can divide the box in slices and calculate the average value of
the order parameter per segment in one slice. The first method gives an
idea of the ordering of a molecule from head to tail, the second method
gives an idea of the ordering as function of the box length.

The electrostatic potential ($\psi$) across the interface can be
computed from a trajectory by evaluating the double integral of the
charge density ($\rho(z)$):

$$
\psi(z) - \psi(-\infty) = - \int_{-\infty}^z dz' \int_{-\infty}^{z'} \rho(z'')dz''/ \epsilon_0
$$

where the position $z=-\infty$ is far enough in the bulk phase such
that the field is zero. With this method, it is possible to “split” the
total potential into separate contributions from lipid and water
molecules. The program `gmx potential ` divides the box
in slices and sums all charges of the atoms in each slice. It then
integrates this charge density to give the electric field, which is in
turn integrated to give the potential. Charge density, electric field,
and potential are written to xvgr input files.

The program `gmx traj ` is a very simple analysis program. All
it does is print the coordinates, velocities, or forces of selected
atoms. It can also calculate the center of mass of one or more molecules
and print the coordinates of the center of mass to three files. By
itself, this is probably not a very useful analysis, but having the
coordinates of selected molecules or atoms can be very handy for further
analysis, not only in interfacial systems.

The program `gmx density ` calculates the mass density of
groups and gives a plot of the density against a box axis. This is
useful for looking at the distribution of groups or atoms across the
interface.
