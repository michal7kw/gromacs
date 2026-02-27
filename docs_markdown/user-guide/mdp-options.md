> **TODO:** Make more cross-references.


# Molecular dynamics parameters (.mdp options)

## General information

Default values are given in parentheses, or listed first among choices.
The first option in the list is always the default option. Units are
given in square brackets. The difference between a dash and an
underscore is ignored.

A `sample mdp file <mdp>` is available. This should be appropriate to
start a normal simulation. Edit it to suit your specific needs and
desires.

### Preprocessing

include

directories to include in your topology. Format:
`-I/home/john/mylib -I../otherlib`


define

defines to pass to the preprocessor, default is no defines. You can use
any defines to control options in your customized topology files.
Options that act on existing `top` file mechanisms include

> `-DFLEXIBLE` will use flexible water instead of rigid water into your
> topology, this can be useful for normal mode analysis.

> `-DPOSRES` will trigger the inclusion of `posre.itp` into your
> topology, used for implementing position restraints.


### Run control

integrator

(Despite the name, this list includes algorithms that are not actually
integrators over time. `integrator=steep` and all entries following it
are in this category)

md

A leap-frog algorithm for integrating Newton's equations of motion.


md-vv

A velocity Verlet algorithm for integrating Newton's equations of
motion. For constant NVE simulations started from corresponding points
in the same trajectory, the trajectories are analytically, but not
binary, identical to the `integrator=md` leap-frog integrator. The
kinetic energy is determined from the whole step velocities and is
therefore slightly too high. The advantage of this integrator is more
accurate, reversible Nose-Hoover and Parrinello-Rahman coupling
integration based on Trotter expansion, as well as (slightly too small)
full step velocity output. This all comes at the cost of extra
computation, especially with constraints and extra communication in
parallel. Note that for nearly all production simulations the
`integrator=md` integrator is accurate enough.


md-vv-avek

A velocity Verlet algorithm identical to `integrator=md-vv`, except that
the kinetic energy is determined as the average of the two half step
kinetic energies as in the `integrator=md` integrator, and this thus
more accurate. With Nose-Hoover and/or Parrinello-Rahman coupling this
comes with a slight increase in computational cost.


sd

An accurate and efficient leap-frog stochastic dynamics integrator. With
constraints, coordinates needs to be constrained twice per integration
step. Depending on the computational cost of the force calculation, this
can take a significant part of the simulation time. The temperature for
one or more groups of atoms (`tc-grps`) is set with `ref-t`, the inverse
friction constant for each group is set with `tau-t`. The parameters
`tcoupl` and `nsttcouple` are ignored. The random generator is
initialized with `ld-seed`. When used as a thermostat, an appropriate
value for `tau-t` is 2 ps, since this results in a friction that is
lower than the internal friction of water, while it is high enough to
remove excess heat NOTE: temperature deviations decay twice as fast as
with a Berendsen thermostat with the same `tau-t`.


bd

An Euler integrator for Brownian or position Langevin dynamics. The
velocity is the force divided by a friction coefficient (`bd-fric`) plus
random thermal noise (`ref-t`). When `bd-fric` is 0, the friction
coefficient for each particle is calculated as mass/ `tau-t`, as for the
integrator `integrator=sd`. The random generator is initialized with
`ld-seed`.


steep

A steepest descent algorithm for energy minimization. The maximum step
size is `emstep`, the tolerance is `emtol`.


cg

A conjugate gradient algorithm for energy minimization, the tolerance is
`emtol`. CG is more efficient when a steepest descent step is done every
once in a while, this is determined by `nstcgsteep`. For a minimization
prior to a normal mode analysis, which requires a very high accuracy,
GROMACS should be compiled in double
precision.


l-bfgs

A quasi-Newtonian algorithm for energy minimization according to the
low-memory Broyden-Fletcher-Goldfarb-Shanno approach. In practice this
seems to converge faster than Conjugate Gradients, but due to the
correction steps necessary it is not (yet) parallelized.


nm

Normal mode analysis is performed on the structure in the `tpr` file.
GROMACS should be compiled in double
precision.


tpi

Test particle insertion. The last molecule in the topology is the test
particle. A trajectory must be provided to `mdrun -rerun`. This
trajectory should not contain the molecule to be inserted. Insertions
are performed `nsteps` times in each frame at random locations and with
random orientations of the molecule. When `nstlist` is larger than one,
`nstlist` insertions are performed in a sphere with radius `rtpi` around
a the same random location using the same pair list. Since pair list
construction is expensive, one can perform several extra insertions with
the same list almost for free. The random seed is set with `ld-seed`.
The temperature for the Boltzmann weighting is set with `ref-t`, this
should match the temperature of the simulation of the original
trajectory. Dispersion correction is implemented correctly for TPI. All
relevant quantities are written to the file specified with `mdrun -tpi`.
The distribution of insertion energies is written to the file specified
with `mdrun -tpid`. No trajectory or energy file is written. Parallel
TPI gives identical results to single-node TPI. For charged molecules,
using PME with a fine grid is most accurate and also efficient, since
the potential in the system only needs to be calculated once per frame.


tpic

Test particle insertion into a predefined cavity location. The procedure
is the same as for `integrator=tpi`, except that one coordinate extra is
read from the trajectory, which is used as the insertion location. The
molecule to be inserted should be centered at 0,0,0.
GROMACS does not do this for you, since for
different situations a different way of centering might be optimal. Also
`rtpi` sets the radius for the sphere around this location. Neighbor
searching is done only once per frame, `nstlist` is not used. Parallel
`integrator=tpic` gives identical results to single-rank
`integrator=tpic`.


mimic

Enable MiMiC QM/MM coupling to run hybrid molecular dynamics. Keey in
mind that its required to launch CPMD compiled with MiMiC as well. In
this mode all options regarding integration (T-coupling, P-coupling,
timestep and number of steps) are ignored as CPMD will do the
integration instead. Options related to forces computation (cutoffs, PME
parameters, etc.) are working as usual. Atom selection to define QM
atoms is read from `QMMM-grps`


tinit

\(0\) \[ps\] starting time for your run (only makes sense for time-based
integrators)


dt

(0.001) \[ps\] time step for integration (only makes sense for
time-based integrators)


nsteps

\(0\) maximum number of steps to integrate or minimize, -1 is no maximum


init-step

\(0\) The starting step. The time at step i in a run is calculated as: t
= `tinit` + `dt` \* (`init-step` + i). The free-energy lambda is
calculated as: lambda = `init-lambda` + `delta-lambda` \* (`init-step` +
i). Also non-equilibrium MD parameters can depend on the step number.
Thus for exact restarts or redoing part of a run it might be necessary
to set `init-step` to the step number of the restart frame.
`gmx convert-tpr` does this automatically.


simulation-part

\(0\) A simulation can consist of multiple parts, each of which has a
part number. This option specifies what that number will be, which helps
keep track of parts that are logically the same simulation. This option
is generally useful to set only when coping with a crashed simulation
where files were lost.


mts

no

Evaluate all forces at every integration step.


yes

Use a multiple timing-stepping integrator to evaluate some forces, as
specified by `mts-level2-forces` every `mts-level2-factor` integration
steps. All other forces are evaluated at every step. MTS is currently
only supported with `integrator=md`.


mts-levels

\(2\) The number of levels for the multiple time-stepping scheme.
Currently only 2 is supported.


mts-level2-forces

(longrange-nonbonded) A list of one or more force groups that will be
evaluated only every `mts-level2-factor` steps. Supported entries are:
`longrange-nonbonded`, `nonbonded`, `pair`, `dihedral`, `angle`, `pull`
and `awh`. With `pair` the listed pair forces (such as 1-4) are
selected. With `dihedral` all dihedrals are selected, including cmap.
All other forces, including all restraints, are evaluated and integrated
every step. When PME or Ewald is used for electrostatics and/or LJ
interactions, `longrange-nonbonded` can not be omitted here.


mts-level2-factor

\(2\) \[steps\] Interval for computing the forces in level 2 of the
multiple time-stepping scheme


mass-repartition-factor

\(1\) \[\] Scales the masses of the lightest atoms in the system by this
factor to the mass mMin. All atoms with a mass lower than mMin also have
their mass set to that mMin. The mass change is subtracted from the mass
of the atom the light atom is bound to. If there is no bound atom a
warning is generated. If there is more than one atom bound an error is
generated. If the mass of the bound atom would become lower than mMin an
error is generated. For typical atomistic systems only the masses of
hydrogens are scaled. With `constraints=h-bonds`, a factor of 3 will
usually enable a time step of 4 fs.


comm-mode

Linear

Remove center of mass translational velocity


Angular

Remove center of mass translational and rotational velocity


Linear-acceleration-correction

Remove center of mass translational velocity. Correct the center of mass
position assuming linear acceleration over `nstcomm` steps. This is
useful for cases where an acceleration is expected on the center of mass
which is nearly constant over `nstcomm` steps. This can occur for
example when pulling on a group using an absolute reference.


None

No restriction on the center of mass motion


nstcomm

\(100\) \[steps\] interval for center of mass motion removal


comm-grps

group(s) for center of mass motion removal, default is the whole system


### Langevin dynamics

bd-fric

\(0\) \[amu ps^{-1}\] Brownian dynamics friction coefficient.
When `bd-fric` is 0, the friction coefficient for each particle is
calculated as mass/ `tau-t`.


ld-seed

(-1) \[integer\] used to initialize random generator for thermal noise
for stochastic and Brownian dynamics. When `ld-seed` is set to -1, a
pseudo random seed is used. When running BD or SD on multiple
processors, each processor uses a seed equal to `ld-seed` plus the
processor number.


### Energy minimization

emtol

(10.0) \[kJ mol^{-1} nm^{-1}\] the minimization is
converged when the maximum force is smaller than this value


emstep

(0.01) \[nm\] initial step-size


nstcgsteep

\(1000\) \[steps\] interval of performing 1 steepest descent step while
doing conjugate gradient energy minimization.


nbfgscorr

\(10\) Number of correction steps to use for L-BFGS minimization. A
higher number is (at least theoretically) more accurate, but slower.


### Shell Molecular Dynamics

When shells or flexible constraints are present in the system the
positions of the shells and the lengths of the flexible constraints are
optimized at every time step until either the RMS force on the shells
and constraints is less than `emtol`, or a maximum number of iterations
`niter` has been reached. Minimization is converged when the maximum
force is smaller than `emtol`. For shell MD this value should be 1.0 at
most.

niter

\(20\) maximum number of iterations for optimizing the shell positions
and the flexible constraints.


fcstep

\(0\) \[ps^{2}\] the step size for optimizing the flexible
constraints. Should be chosen as mu/(d2V/dq2) where mu is the reduced
mass of two particles in a flexible constraint and d2V/dq2 is the second
derivative of the potential in the constraint direction. Hopefully this
number does not differ too much between the flexible constraints, as the
number of iterations and thus the runtime is very sensitive to fcstep.
Try several values!


### Test particle insertion

rtpi

(0.05) \[nm\] the test particle insertion radius, see integrators
`integrator=tpi` and `integrator=tpic`


### Output control

nstxout

\(0\) \[steps\] number of steps that elapse between writing coordinates
to the output trajectory file (`trr`), the first and last coordinates
are always written unless 0, which means coordinates are not written
into the trajectory file.


nstvout

\(0\) \[steps\] number of steps that elapse between writing velocities
to the output trajectory file (`trr`), the first and last velocities are
always written unless 0, which means velocities are not written into the
trajectory file.


nstfout

\(0\) \[steps\] number of steps that elapse between writing forces to
the output trajectory file (`trr`), the first and last forces are always
written, unless 0, which means forces are not written into the
trajectory file.


nstlog

\(1000\) \[steps\] number of steps that elapse between writing energies
to the log file, the first and last energies are always written.


nstcalcenergy

\(100\) number of steps that elapse between calculating the energies, 0
is never. This option is only relevant with dynamics. This option
affects the performance in parallel simulations, because calculating
energies requires global communication between all processes which can
become a bottleneck at high parallelization.


nstenergy

\(1000\) \[steps\] number of steps that elapse between writing energies
to the energy file (`edr`), the first and last energies are always
written, should be a multiple of `nstcalcenergy`. Note that the exact
sums and fluctuations over all MD steps modulo `nstcalcenergy` are
stored in the energy file, so `gmx energy` can report exact energy
averages and fluctuations also when `nstenergy` \> 1


nstxout-compressed

\(0\) \[steps\] number of steps that elapse between writing position
coordinates using lossy compression (`xtc` file), the first and last
coordinates are always written, unless 0, which means that there is no
compressed coordinates output.


compressed-x-precision

\(1000\) \[real\] precision with which to write to the compressed
trajectory file


compressed-x-grps

group(s) to write to the compressed trajectory file, by default the
whole system is written (if `nstxout-compressed` \> 0)


energygrps

group(s) for which to write to write short-ranged non-bonded potential
energies to the energy file (not supported on GPUs)


### Neighbor searching

cutoff-scheme

Verlet

Generate a pair list with buffering. The buffer size is automatically
set based on `verlet-buffer-tolerance`, unless this is set to -1, in
which case `rlist` will be used.


group

Generate a pair list for groups of atoms, corresponding to the charge
groups in the topology. This option is no longer supported.


nstlist

10) \[steps\]

\>0

Interval between steps that update the neighbor list. When dynamics and
`verlet-buffer-tolerance` set, `nstlist` is actually a minimum value and
`gmx mdrun` might increase it, unless it is set to 1. With parallel
simulations and/or non-bonded force calculation on the GPU, a value of
20 or 40 often gives the best performance. With energy minimization this
parameter is not used as the pair list is updated when at least one atom
has moved by more than half the pair list buffer size.


0

The neighbor list is only constructed once and never updated. This is
mainly useful for vacuum simulations in which all particles see each
other. But vacuum simulations are (temporarily) not supported.


\<0

Unused.


pbc

xyz

Use periodic boundary conditions in all directions.


no

Use no periodic boundary conditions, ignore the box. To simulate without
cut-offs, set all cut-offs and `nstlist` to 0. For best performance
without cut-offs on a single MPI rank, set `nstlist` to zero.


xy

Use periodic boundary conditions in x and y directions only. This can be
used in combination with [walls](#walls). Without walls or with only one
wall the system size is infinite in the z direction. Therefore pressure
coupling or Ewald summation methods can not be used. These disadvantages
do not apply when two walls are used.


periodic-molecules

no

molecules are finite, fast molecular PBC can be used


yes

for systems with molecules that couple to themselves through the
periodic boundary conditions, this requires a slower PBC algorithm and
molecules are not made whole in the output


verlet-buffer-tolerance

(0.005) \[kJ mol^{-1} ps^{-1}\]

Used when performing a simulation with dynamics. This sets the maximum
allowed error for pair interactions per particle caused by the Verlet
buffer, which indirectly sets `rlist`. As both `nstlist` and the Verlet
buffer size are fixed (for performance reasons), particle pairs not in
the pair list can occasionally get within the cut-off distance during
`nstlist` -1 steps. This causes very small jumps in the energy. In a
constant-temperature ensemble, these very small energy jumps can be
estimated for a given cut-off and `rlist`. The estimate assumes a
homogeneous particle distribution, hence the errors might be slightly
underestimated for multi-phase systems. (See the [reference manual]()
for details). For longer pair-list life-time (`nstlist` -1) \* `dt` the
buffer is overestimated, because the interactions between particles are
ignored. Combined with cancellation of errors, the actual drift of the
total energy is usually one to two orders of magnitude smaller. Note
that the generated buffer size takes into account that the
GROMACS pair-list setup leads to a reduction
in the drift by a factor 10, compared to a simple particle-pair based
list. Without dynamics (energy minimization etc.), the buffer is 5% of
the cut-off. For NVE simulations the initial temperature is used, unless
this is zero, in which case a buffer of 10% is used. For NVE simulations
the tolerance usually needs to be lowered to achieve proper energy
conservation on the nanosecond time scale. To override the automated
buffer setting, use `verlet-buffer-tolerance` =-1 and set `rlist`
manually.


verlet-buffer-pressure-tolerance

(0.5) \[bar\]

Used when performing a simulation with dynamics and only active when
`verlet-buffer-tolerance` is positive. This sets the maximum tolerated
error in the average pressure due to missing Lennard-Jones interactions
of particle pairs that are not in the pair list, but come within `rvdw`
range as the pair list ages. As for the drift tolerance, the
(over)estimate of the pressure error is tight at short times. At longer
time it turns into a significant overestimate, because interactions
limit the displacement of particles. Note that the default tolerance of
0.5 bar corresponds to a maximum relative deviation of the density of
liquid water of 2e-5.


rlist

\(1\) \[nm\] Cut-off distance for the short-range neighbor list. With
dynamics, this is by default set by the `verlet-buffer-tolerance` and
`verlet-buffer-pressure-tolerance` options and the value of `rlist` is
ignored. Without dynamics, this is by default set to the maximum cut-off
plus 5% buffer, except for test particle insertion, where the buffer is
managed exactly and automatically. For NVE simulations, where the
automated setting is not possible, the advised procedure is to run
`gmx grompp` with an NVT setup with the expected temperature and copy
the resulting value of `rlist` to the NVE setup.


### Electrostatics

coulombtype

Cut-off

Plain cut-off with pair list radius `rlist` and Coulomb cut-off
`rcoulomb`, where `rlist` \>= `rcoulomb`. Note that with the (default)
setting of `coulomb-modifier` =Potential-shift not only the potentials
between interacting pairs are shifted to be zero at the cut-off, but the
same shift is also applied to excluded pairs. This does not lead to
forces between excluded pairs, but does add a constant offset to the
total Coulomb potential.


Ewald

Classical Ewald sum electrostatics. The real-space cut-off `rcoulomb`
should be equal to `rlist`. Use *e.g.* `rlist` =0.9, `rcoulomb` =0.9.
The highest magnitude of wave vectors used in reciprocal space is
controlled by `fourierspacing`. The relative accuracy of
direct/reciprocal space is controlled by `ewald-rtol`.

NOTE: Ewald scales as O(N^{3/2}) and is thus extremely slow for
large systems. It is included mainly for reference - in most cases PME
will perform much better.


PME

Fast smooth Particle-Mesh Ewald (SPME) electrostatics. Direct space is
similar to the Ewald sum, while the reciprocal part is performed with
FFTs. Grid dimensions are controlled with `fourierspacing` and the
interpolation order with `pme-order`. With a grid spacing of 0.1 nm and
cubic interpolation the electrostatic forces have an accuracy of
2-3\*10^{-4}. Since the error from the vdw-cutoff is larger than
this you might try 0.15 nm. When running in parallel the interpolation
parallelizes better than the FFT, so try decreasing grid dimensions
while increasing interpolation.


P3M-AD

Particle-Particle Particle-Mesh algorithm with analytical derivative for
for long-range electrostatic interactions. The method and code is
identical to SPME, except that the influence function is optimized for
the grid. This gives a slight increase in accuracy.


Reaction-Field

Reaction field electrostatics with Coulomb cut-off `rcoulomb`, where
`rlist` \>= `rvdw`. The dielectric constant beyond the cut-off is
`epsilon-rf`. The dielectric constant can be set to infinity by setting
`epsilon-rf` =0.


User

Currently unsupported. `gmx mdrun` will now expect to find a file
`table.xvg` with user-defined potential functions for repulsion,
dispersion and Coulomb. When pair interactions are present, `gmx
mdrun` also expects to find a file `tablep.xvg` for the pair
interactions. When the same interactions should be used for non-bonded
and pair interactions the user can specify the same file name for both
table files. These files should contain 7 columns: the `x` value,
`f(x)`, `-f'(x)`, `g(x)`, `-g'(x)`, `h(x)`, `-h'(x)`, where `f(x)` is
the Coulomb function, `g(x)` the dispersion function and `h(x)` the
repulsion function. When `vdwtype` is not set to User the values for
`g`, `-g'`, `h` and `-h'` are ignored. For the non-bonded interactions
`x` values should run from 0 to the largest cut-off distance +
`table-extension` and should be uniformly spaced. For the pair
interactions the table length in the file will be used. The optimal
spacing, which is used for non-user tables, is `0.002 nm` when you run
in mixed precision or `0.0005 nm` when you run in double precision. The
function value at `x=0` is not important. More information is in the
printed manual.


PME-Switch

Currently unsupported. A combination of PME and a switch function for
the direct-space part (see above). `rcoulomb` is allowed to be smaller
than `rlist`.


PME-User

Currently unsupported. A combination of PME and user tables (see above).
`rcoulomb` is allowed to be smaller than `rlist`. The PME mesh
contribution is subtracted from the user table by `gmx mdrun`. Because
of this subtraction the user tables should contain about 10 decimal
places.


PME-User-Switch

Currently unsupported. A combination of PME-User and a switching
function (see above). The switching function is applied to final
particle-particle interaction, *i.e.* both to the user supplied function
and the PME Mesh correction part.


coulomb-modifier

Potential-shift

Shift the Coulomb potential by a constant such that it is zero at the
cut-off. This makes the potential the integral of the force. Note that
this does not affect the forces or the sampling.


None

Use an unmodified Coulomb potential. This can be useful when comparing
energies with those computed with other software.


rcoulomb-switch

\(0\) \[nm\] where to start switching the Coulomb potential, only
relevant when force or potential switching is used


rcoulomb

\(1\) \[nm\] The distance for the Coulomb cut-off. Note that with PME
this value can be increased by the PME tuning in `gmx mdrun` along with
the PME grid spacing.


epsilon-r

\(1\) The relative dielectric constant. A value of 0 means infinity.


epsilon-rf

\(0\) The relative dielectric constant of the reaction field. This is
only used with reaction-field electrostatics. A value of 0 means
infinity.


### Van der Waals

vdwtype

Cut-off

Plain cut-off with pair list radius `rlist` and VdW cut-off `rvdw`,
where `rlist` \>= `rvdw`.


PME

Fast smooth Particle-mesh Ewald (SPME) for VdW interactions. The grid
dimensions are controlled with `fourierspacing` in the same way as for
electrostatics, and the interpolation order is controlled with
`pme-order`. The relative accuracy of direct/reciprocal space is
controlled by `ewald-rtol-lj`, and the specific combination rules that
are to be used by the reciprocal routine are set using
`lj-pme-comb-rule`.


Shift

This functionality is deprecated and replaced by using `vdwtype=Cut-off`
with `vdw-modifier=Force-switch`. The LJ (not Buckingham) potential is
decreased over the whole range and the forces decay smoothly to zero
between `rvdw-switch` and `rvdw`.


Switch

This functionality is deprecated and replaced by using `vdwtype=Cut-off`
with `vdw-modifier=Potential-switch`. The LJ (not Buckingham) potential
is normal out to `rvdw-switch`, after which it is switched off to reach
zero at `rvdw`. Both the potential and force functions are continuously
smooth, but be aware that all switch functions will give rise to a bulge
(increase) in the force (since we are switching the potential).


User

Currently unsupported. See `coulombtype=User` for instructions. The
function value at zero is not important. When you want to use LJ
correction, make sure that `rvdw` corresponds to the cut-off in the
user-defined function. When `coulombtype` is not set to User the values
for the `f` and `-f'` columns are ignored.


vdw-modifier

Potential-shift

Shift the Van der Waals potential by a constant such that it is zero at
the cut-off. This makes the potential the integral of the force. Note
that this does not affect the forces or the sampling.


None

Use an unmodified Van der Waals potential. This can be useful when
comparing energies with those computed with other software.


Force-switch

Smoothly switches the forces to zero between `rvdw-switch` and `rvdw`.
This shifts the potential shift over the whole range and switches it to
zero at the cut-off. Note that this is more expensive to calculate than
a plain cut-off and it is not required for energy conservation, since
Potential-shift conserves energy just as well.


Potential-switch

Smoothly switches the potential to zero between `rvdw-switch` and
`rvdw`. Note that this introduces articifically large forces in the
switching region and is much more expensive to calculate. This option
should only be used if the force field you are using requires this.


rvdw-switch

\(0\) \[nm\] where to start switching the LJ force and possibly the
potential, only relevant when force or potential switching is used


rvdw

\(1\) \[nm\] distance for the LJ or Buckingham cut-off


DispCorr

no

do not apply any correction


EnerPres

apply long-range dispersion corrections for Energy and Pressure


Ener

apply long-range dispersion corrections for Energy only


### Tables

table-extension

\(1\) \[nm\] Extension of the non-bonded potential lookup tables beyond
the largest cut-off distance. With actual non-bonded interactions the
tables are never accessed beyond the cut-off. But a longer table length
might be needed for the 1-4 interactions, which are always tabulated
irrespective of the use of tables for the non-bonded interactions.


energygrp-table

Currently unsupported. When user tables are used for electrostatics
and/or VdW, here one can give pairs of energy groups for which separate
user tables should be used. The two energy groups will be appended to
the table file name, in order of their definition in `energygrps`,
separated by underscores. For example, if `energygrps = Na Cl Sol` and
`energygrp-table = Na Na Na Cl`, `gmx mdrun` will read `table_Na_Na.xvg`
and `table_Na_Cl.xvg` in addition to the normal `table.xvg` which will
be used for all other energy group pairs.


### Ewald

fourierspacing

(0.12) \[nm\] For ordinary Ewald, the ratio of the box dimensions and
the spacing determines a lower bound for the number of wave vectors to
use in each (signed) direction. For PME and P3M, that ratio determines a
lower bound for the number of Fourier-space grid points that will be
used along that axis. In all cases, the number for each direction can be
overridden by entering a non-zero value for that `fourier-nx` direction.
For optimizing the relative load of the particle-particle interactions
and the mesh part of PME, it is useful to know that the accuracy of the
electrostatics remains nearly constant when the Coulomb cut-off and the
PME grid spacing are scaled by the same factor. Note that this spacing
can be scaled up along with `rcoulomb` by the PME tuning in `gmx mdrun`.


fourier-nx


fourier-ny


fourier-nz

\(0\) Highest magnitude of wave vectors in reciprocal space when using
Ewald. Grid size when using PME or P3M. These values override
`fourierspacing` per direction. The best choice is powers of 2, 3, 5 and
7. Avoid large primes. Note that these grid sizes can be reduced along
with scaling up `rcoulomb` by the PME tuning in `gmx mdrun`.


pme-order

\(4\) The number of grid points along a dimension to which a charge is
mapped. The actual order of the PME interpolation is one less, e.g. the
default of 4 gives cubic interpolation. Supported values are 3 to 12
(max 8 for P3M-AD). When running in parallel, it can be worth to switch
to 5 and simultaneously increase the grid spacing. Note that on the CPU
only values 4 and 5 have SIMD acceleration and GPUs only support the
value 4.


ewald-rtol

(10^{-5}) The relative strength of the Ewald-shifted direct
potential at `rcoulomb` is given by `ewald-rtol`. Decreasing this will
give a more accurate direct sum, but then you need more wave vectors for
the reciprocal sum.


ewald-rtol-lj

(10^{-3}) When doing PME for VdW-interactions, `ewald-rtol-lj`
is used to control the relative strength of the dispersion potential at
`rvdw` in the same way as `ewald-rtol` controls the electrostatic
potential.


lj-pme-comb-rule

(Geometric) The combination rules used to combine VdW-parameters in the
reciprocal part of LJ-PME. Geometric rules are much faster than
Lorentz-Berthelot and usually the recommended choice, even when the rest
of the force field uses the Lorentz-Berthelot rules.

Geometric

Apply geometric combination rules


Lorentz-Berthelot

Apply Lorentz-Berthelot combination rules


ewald-geometry

3d

The Ewald sum is performed in all three dimensions.


3dc

The reciprocal sum is still performed in 3D, but a force and potential
correction applied in the `z` dimension to produce a pseudo-2D
summation. If your system has a slab geometry in the `x-y` plane you can
try to increase the `z`-dimension of the box (a box height of 3 times
the slab height is usually ok) and use this option.


epsilon-surface

\(0\) This controls the dipole correction to the Ewald summation in 3D.
The default value of zero means it is turned off. Turn it on by setting
it to the value of the relative permittivity of the imaginary surface
around your infinite system. Be careful - you should not use this if you
have free mobile charges in your system. This value does not affect the
slab 3DC variant of the long-range corrections.


### Temperature coupling

ensemble-temperature-setting

auto

With this setting `gmx grompp` will determine which of the next three
settings is available and choose the appropriate one. When all atoms are
coupled to a temperature bath with the same temperature, a constant
ensemble temperature is chosen and the value is taken from the
temperature bath.


constant

The system has a constant ensemble temperature given by
`ensemble-temperature`. A constant ensemble temperature is required for
certain sampling algorithms such as AWH.


variable

The system has a variable ensemble temperature due to simulated
annealing or simulated tempering. The system ensemble temperature is set
dynamically during the simulation.


not-available

The system has no ensemble temperature.


ensemble-temperature

(-1) \[K\]

The ensemble temperature for the system. The input value is only used
with `ensemble-temperature-setting=constant`. By default the ensemble
temperature is copied from the temperature of the thermal bath (when
used).


tcoupl

no

No temperature coupling.


berendsen

Temperature coupling with a Berendsen thermostat to a bath with
temperature `ref-t`, with time constant `tau-t`. Several groups can be
coupled separately, these are specified in the `tc-grps` field separated
by spaces. This is a historical thermostat needed to be able to
reproduce previous simulations, but we strongly recommend not to use it
for new production runs. Consult the manual for details.


nose-hoover

Temperature coupling using a Nose-Hoover extended ensemble. The
reference temperature and coupling groups are selected as above, but in
this case `tau-t` controls the period of the temperature fluctuations at
equilibrium, which is slightly different from a relaxation time. For NVT
simulations the conserved energy quantity is written to the energy and
log files.


andersen

Temperature coupling by randomizing a fraction of the particle
velocities at each timestep. Reference temperature and coupling groups
are selected as above. `tau-t` is the average time between randomization
of each molecule. Inhibits particle dynamics somewhat, but has little or
no ergodicity issues. Currently only implemented with velocity Verlet,
and not implemented with constraints.


andersen-massive

Temperature coupling by randomizing velocities of all particles at
infrequent timesteps. Reference temperature and coupling groups are
selected as above. `tau-t` is the time between randomization of all
molecules. Inhibits particle dynamics somewhat, but has little or no
ergodicity issues. Currently only implemented with velocity Verlet.


v-rescale

Temperature coupling using velocity rescaling with a stochastic term
(JCP 126, 014101). This thermostat is similar to Berendsen coupling,
with the same scaling using `tau-t`, but the stochastic term ensures
that a proper canonical ensemble is generated. The random seed is set
with `ld-seed`. This thermostat works correctly even for `tau-t` =0. For
NVT simulations the conserved energy quantity is written to the energy
and log file.


nsttcouple

(-1) The interval between steps that couple the temperature. The default
value of -1 sets `nsttcouple` equal to 100, or fewer steps if required
for accurate integration (5 steps per tau for first order coupling, 20
steps per tau for second order coupling). Note that the default value is
large in order to reduce the overhead of the additional computation and
communication required for obtaining the kinetic energy. For velocity
Verlet integrators `nsttcouple` is set to 1.


nh-chain-length

\(10\) The number of chained Nose-Hoover thermostats for velocity Verlet
integrators, the leap-frog `integrator=md` integrator only supports 1.
Data for the NH chain variables is not printed to the `edr` file by
default, but can be turned on with the
`print-nose-hoover-chain-variables` option.


print-nose-hoover-chain-variables

no

Do not store Nose-Hoover chain variables in the energy file.


yes

Store all positions and velocities of the Nose-Hoover chain in the
energy file.


tc-grps

groups to couple to separate temperature baths


tau-t

\[ps\] time constant for coupling (one for each group in `tc-grps`), -1
means no temperature coupling


ref-t

\[K\] reference temperature for coupling (one for each group in
`tc-grps`)


### Pressure coupling

pcoupl

no

No pressure coupling. This means a fixed box size.


Berendsen

Exponential relaxation pressure coupling with time constant `tau-p`. The
box is scaled every `nstpcouple` steps. This barostat does not yield a
correct thermodynamic ensemble; it is only included to be able to
reproduce previous runs, and we strongly recommend against using it for
new simulations. See the manual for details.


C-rescale

Exponential relaxation pressure coupling with time constant `tau-p`,
including a stochastic term to enforce correct volume fluctuations. The
box is scaled every `nstpcouple` steps. It can be used for both
equilibration and production.


Parrinello-Rahman

Extended-ensemble pressure coupling where the box vectors are subject to
an equation of motion. The equation of motion for the atoms is coupled
to this. No instantaneous scaling takes place. As for Nose-Hoover
temperature coupling the time constant `tau-p` is the period of pressure
fluctuations at equilibrium. This is a good method when you want to
apply pressure scaling during data collection, but beware that you can
get very large oscillations if you are starting from a different
pressure. For simulations where the exact fluctations of the NPT
ensemble are important, or if the pressure coupling time is very short,
it may not be appropriate, as the previous time step pressure is used in
some steps of the GROMACS implementation for
the current time step pressure.


MTTK

Martyna-Tuckerman-Tobias-Klein implementation, only useable with
`integrator=md-vv` or `integrator=md-vv-avek`, very similar to
Parrinello-Rahman. As for Nose-Hoover temperature coupling the time
constant `tau-p` is the period of pressure fluctuations at equilibrium.
This is probably a better method when you want to apply pressure scaling
during data collection, but beware that you can get very large
oscillations if you are starting from a different pressure. This
requires a constant ensemble temperature for the system. It only
supports isotropic scaling, and only works without constraints. MTTK
coupling is deprecated.


pcoupltype

Specifies the kind of isotropy of the pressure coupling used. Each kind
takes one or more values for `compressibility` and `ref-p`. Only a
single value is permitted for `tau-p`.

isotropic

Isotropic pressure coupling with time constant `tau-p`. One value each
for `compressibility` and `ref-p` is required.


semiisotropic

Pressure coupling which is isotropic in the `x` and `y` direction, but
different in the `z` direction. This can be useful for membrane
simulations. Two values each for `compressibility` and `ref-p` are
required, for `x/y` and `z` directions respectively.


anisotropic

Same as before, but 6 values are needed for `xx`, `yy`, `zz`, `xy/yx`,
`xz/zx` and `yz/zy` components, respectively. When the off-diagonal
compressibilities are set to zero, a rectangular box will stay
rectangular. Beware that anisotropic scaling can lead to extreme
deformation of the simulation box.


surface-tension

Surface tension coupling for surfaces parallel to the xy-plane. Uses
normal pressure coupling for the `z`-direction, while the surface
tension is coupled to the `x/y` dimensions of the box. The first `ref-p`
value is the reference surface tension times the number of surfaces
`bar nm`, the second value is the reference `z`-pressure `bar`. The two
`compressibility` values are the compressibility in the `x/y` and `z`
direction respectively. The value for the `z`-compressibility should be
reasonably accurate since it influences the convergence of the
surface-tension, it can also be set to zero to have a box with constant
height.


nstpcouple

(-1) The interval between steps that couple the pressure. The default
value of -1 sets `nstpcouple` equal to 100, or fewer steps if required
for accurate integration (5 steps per tau for first order coupling, 20
steps per tau for second order coupling). Note that the default value is
large in order to reduce the overhead of the additional computation and
communication required for obtaining the virial and kinetic energy. For
velocity Verlet integrators `nsttcouple` is set to 1.


tau-p

\(5\) \[ps\] The time constant for pressure coupling (one value for all
directions).


compressibility

\[bar^{-1}\] The compressibility (NOTE: this is now really in
bar^{-1}) For water at 1 atm and 300 K the compressibility is
4.5e-5 bar^{-1}. The number of required values is implied by
`pcoupltype`.


ref-p

\[bar\] The reference pressure for coupling. The number of required
values is implied by `pcoupltype`.


refcoord-scaling

no

The reference coordinates for position restraints are not modified. Note
that with this option the virial and pressure might be ill defined, see
`here <reference-manual-position-restraints>` for more details.


all

The reference coordinates are scaled with the scaling matrix of the
pressure coupling.


com

Scale the center of mass of the reference coordinates with the scaling
matrix of the pressure coupling. The vectors of each reference
coordinate to the center of mass are not scaled. Only one COM is used,
even when there are multiple molecules with position restraints. For
calculating the COM of the reference coordinates in the starting
configuration, periodic boundary conditions are not taken into account.
Note that with this option the virial and pressure might be ill defined,
see `here <reference-manual-position-restraints>` for more details.


### Simulated annealing

Simulated annealing is controlled separately for each temperature group
in GROMACS. The reference temperature is a
piecewise linear function, but you can use an arbitrary number of points
for each group, and choose either a single sequence or a periodic
behaviour for each group. The actual annealing is performed by
dynamically changing the reference temperature used in the thermostat
algorithm selected, so remember that the system will usually not
instantaneously reach the reference temperature!

annealing

Type of annealing for each temperature group

no

No simulated annealing - just couple to reference temperature value.


single

A single sequence of annealing points. If your simulation is longer than
the time of the last point, the temperature will be coupled to this
constant value after the annealing sequence has reached the last time
point.


periodic

The annealing will start over at the first reference point once the last
reference time is reached. This is repeated until the simulation ends.


annealing-npoints

A list with the number of annealing reference/control points used for
each temperature group. Use 0 for groups that are not annealed. The
number of entries should equal the number of temperature groups.


annealing-time

List of times at the annealing reference/control points for each group.
If you are using periodic annealing, the times will be used modulo the
last value, *i.e.* if the values are 0, 5, 10, and 15, the coupling will
restart at the 0ps value after 15ps, 30ps, 45ps, etc. The number of
entries should equal the sum of the numbers given in
`annealing-npoints`.


annealing-temp

List of temperatures at the annealing reference/control points for each
group. The number of entries should equal the sum of the numbers given
in `annealing-npoints`.


Confused? OK, let's use an example. Assume you have two temperature
groups, set the group selections to `annealing = single periodic`, the
number of points of each group to `annealing-npoints = 3 4`, the times
to `annealing-time = 0 3 6 0 2 4 6` and finally temperatures to
`annealing-temp = 298 280 270 298 320 320 298`. The first group will be
coupled to 298K at 0ps, but the reference temperature will drop linearly
to reach 280K at 3ps, and then linearly between 280K and 270K from 3ps
to 6ps. After this is stays constant, at 270K. The second group is
coupled to 298K at 0ps, it increases linearly to 320K at 2ps, where it
stays constant until 4ps. Between 4ps and 6ps it decreases to 298K, and
then it starts over with the same pattern again, *i.e.* rising linearly
from 298K to 320K between 6ps and 8ps. Check the summary printed by
`gmx grompp` if you are unsure!

### Velocity generation

gen-vel

no

Do not generate velocities. The velocities are set to zero when there
are no velocities in the input structure file.


yes

Generate velocities in `gmx grompp` according to a Maxwell distribution
at temperature `gen-temp`, with random seed `gen-seed`. This is only
meaningful with `integrator=md`.


gen-temp

\(300\) \[K\] temperature for Maxwell distribution


gen-seed

(-1) \[integer\] used to initialize random generator for random
velocities, when `gen-seed` is set to -1, a pseudo random seed is used.


### Bonds

constraints

Controls which bonds in the topology will be converted to rigid
holonomic constraints. Note that typical rigid water models do not have
bonds, but rather a specialized `[settles]` directive, so are not
affected by this keyword.

none

No bonds converted to constraints.


h-bonds

Convert the bonds with H-atoms to constraints.


all-bonds

Convert all bonds to constraints. Note that many force fields have been
parameterized with flexible bonds between heavy atoms and that
constraining these bonds can increase energy barriers of dihedrals.


h-angles

Convert all bonds to constraints and convert the angles that involve
H-atoms to bond-constraints. This option is deprecated.


all-angles

Convert all bonds to constraints and all angles to bond-constraints.
This option is deprecated.


constraint-algorithm

Chooses which solver satisfies any non-SETTLE holonomic constraints.

LINCS

LINear Constraint Solver. With domain decomposition the parallel version
P-LINCS is used. The accuracy in set with `lincs-order`, which sets the
number of matrices in the expansion for the matrix inversion. After the
matrix inversion correction the algorithm does an iterative correction
to compensate for lengthening due to rotation. The number of such
iterations can be controlled with `lincs-iter`. The root mean square
relative constraint deviation is printed to the log file every `nstlog`
steps. If a bond rotates more than `lincs-warnangle` in one step, a
warning will be printed both to the log file and to `stderr`. LINCS
should not be used with coupled angle constraints.


SHAKE

SHAKE is slightly slower and less stable than LINCS, but does work with
angle constraints. The relative tolerance is set with `shake-tol`,
0.0001 is a good value for "normal" MD. SHAKE does not support
constraints between atoms on different decomposition domains, so it can
only be used with domain decomposition when so-called update-groups are
used, which is usually the case when only bonds involving hydrogens are
constrained. SHAKE can not be used with energy minimization.


continuation

This option was formerly known as `unconstrained-start`.

no

apply constraints to the start configuration and reset shells


yes

do not apply constraints to the start configuration and do not reset
shells, useful for exact continuation and reruns


shake-tol

(0.0001) relative tolerance for SHAKE


lincs-order

\(4\) Highest order in the expansion of the constraint coupling matrix.
When constraints form triangles, an additional expansion of the same
order is applied on top of the normal expansion only for the couplings
within such triangles. For "normal" MD simulations an order of 4 usually
suffices, 6 is needed for large time-steps with virtual sites or BD. For
accurate energy minimization in double precision an order of 8 or more
might be required. Note that in single precision an order higher than 6
will often lead to worse accuracy due to amplification of rounding
errors. With domain decomposition, the cell size is limited by the
distance spanned by `lincs-order` +1 constraints. When one wants to
scale further than this limit, one can decrease `lincs-order` and
increase `lincs-iter`, since the accuracy does not deteriorate when (1+
`lincs-iter` )\* `lincs-order` remains constant.


lincs-iter

\(1\) Number of iterations to correct for rotational lengthening in
LINCS. For normal runs a single step is sufficient, but for NVE runs
where you want to conserve energy accurately or for accurate energy
minimization in double precision you might want to increase it to 2.
Note that in single precision using more than 1 iteration will often
lead to worse accuracy due to amplification of rounding errors.


lincs-warnangle

\(30\) \[deg\] maximum angle that a bond can rotate before LINCS will
complain


morse

no

bonds are represented by a harmonic potential


yes

bonds are represented by a Morse potential


### Energy group exclusions

energygrp-excl

Exclusion between pairs of energy groups are currently not supported.


### Walls

nwall

\(0\) When set to 1 there is a wall at `z=0`, when set to 2 there is
also a wall at `z=z-box`. Walls can only be used with `pbc` `=xy`. When
set to 2, pressure coupling and Ewald summation can be used (it is
usually best to use semiisotropic pressure coupling with the `x/y`
compressibility set to 0, as otherwise the surface area will change).
Walls interact wit the rest of the system through an optional
`wall-atomtype`. Energy groups `wall0` and `wall1` (for `nwall` =2) are
added automatically to monitor the interaction of energy groups with
each wall. The center of mass motion removal will be turned off in the
`z`-direction.


wall-atomtype

the atom type name in the force field for each wall. By (for example)
defining a special wall atom type in the topology with its own
combination rules, this allows for independent tuning of the interaction
of each atomtype with the walls.


wall-type

9-3

LJ integrated over the volume behind the wall: 9-3 potential


10-4

LJ integrated over the wall surface: 10-4 potential


12-6

direct LJ potential with the `z` distance from the wall


table

user-defined potentials indexed with the `z` distance from the wall, the
tables are read analogously to the `energygrp-table` option, where the
first name is for a "normal" energy group and the second name is `wall0`
or `wall1`, only the dispersion and repulsion columns are used


wall-r-linpot

(-1) \[nm\] Below this distance from the wall the potential is continued
linearly and thus the force is constant. Setting this option to a
postive value is especially useful for equilibration when some atoms are
beyond a wall. When the value is \<=0 (\<0 for `wall-type` =table), a
fatal error is generated when atoms are beyond a wall.


wall-density

\[nm^{-3}\] / \[nm^{-2}\] the number density of the
atoms for each wall for wall types 9-3 and 10-4


wall-ewald-zfac

\(3\) The scaling factor for the third box vector for Ewald summation
only, the minimum is 2. Ewald summation can only be used with `nwall`
=2, where one should use `ewald-geometry` `=3dc`. The empty layer in the
box serves to decrease the unphysical Coulomb interaction between
periodic images.


### COM pulling

Sets whether pulling on collective variables is active. Note that where
pulling coordinates are applicable, there can be more than one (set with
`pull-ncoords`) and multiple related `mdp` variables will exist
accordingly. Documentation references to things like `pull-coord1-vec`
should be understood to apply to to the applicable pulling coordinate,
eg. the second pull coordinate is described by pull-coord2-vec,
pull-coord2-k, and so on.

pull

no

No center of mass pulling. All the following pull options will be
ignored (and if present in the `mdp` file, they unfortunately generate
warnings)


yes

Center of mass pulling will be applied on 1 or more groups using 1 or
more pull coordinates.


pull-cylinder-r

(1.5) \[nm\] the radius of the cylinder for
`pull-coord1-geometry=cylinder`


pull-constr-tol

(10^{-6}) the relative constraint tolerance for constraint
pulling


pull-print-com

no

do not print the COM for any group


yes

print the COM of all groups for all pull coordinates to the `pullx.xvg`
file.


pull-print-ref-value

no

do not print the reference value for each pull coordinate


yes

print the reference value for each pull coordinate to the `pullx.xvg`
file.


pull-print-components

no

only print the distance for each pull coordinate


yes

print the distance and Cartesian components selected in
`pull-coord1-dim` to the `pullx.xvg` file.


pull-nstxout

\(50\) interval for writing out the COMs of all the pull groups (0 is
never) to the `pullx.xvg` file.


pull-nstfout

\(50\) interval for writing out the force of all the pulled groups (0 is
never) to the `pullf.xvg` file.


pull-pbc-ref-prev-step-com

no

Use the reference atom (`pull-group1-pbcatom`) for the treatment of
periodic boundary conditions.


yes

Use the COM of the previous step as reference for the treatment of
periodic boundary conditions. The reference is initialized using the
reference atom (`pull-group1-pbcatom`), which should be located
centrally in the group. Using the COM from the previous step can be
useful if one or more pull groups are large or very flexible.


pull-xout-average

no

Write the instantaneous coordinates for all the pulled groups.


yes

Write the average coordinates (since last output) for all the pulled
groups. N.b., some analysis tools might expect instantaneous pull
output.


pull-fout-average

no

Write the instantaneous force for all the pulled groups.


yes

Write the average force (since last output) for all the pulled groups.
N.b., some analysis tools might expect instantaneous pull output.


pull-ngroups

\(1\) The number of pull groups, not including the absolute reference
group, when used. Pull groups can be reused in multiple pull
coordinates. Below only the pull options for group 1 are given, further
groups simply increase the group index number.


pull-ncoords

\(1\) The number of pull coordinates. Below only the pull options for
coordinate 1 are given, further coordinates simply increase the
coordinate index number.


pull-group1-name

The name of the pull group, is looked up in the index file or in the
default groups to obtain the atoms involved.


pull-group1-weights

Optional relative weights which are multiplied with the masses of the
atoms to give the total weight for the COM. The number of weights should
be 0, meaning all 1, or the number of atoms in the pull group.


pull-group1-pbcatom

\(0\) The reference atom for the treatment of periodic boundary
conditions inside the group (this has no effect on the treatment of the
pbc between groups). This option is only important when the diameter of
the pull group is larger than half the shortest box vector. For
determining the COM, all atoms in the group are put at their periodic
image which is closest to `pull-group1-pbcatom`. A value of 0 means that
the middle atom (number wise) is used, which is only safe for small
groups. `gmx grompp` checks that the maximum distance from the reference
atom (specifically chosen, or not) to the other atoms in the group is
not too large. This parameter is not used with `pull-coord1-geometry`
cylinder. A value of -1 turns on cosine weighting, which is useful for a
group of molecules in a periodic system, *e.g.* a water slab (see Engin
et al. J. Chem. Phys. B 2010).


pull-coord1-type

umbrella

Center of mass pulling using an umbrella potential between the reference
group and one or more groups.


constraint

Center of mass pulling using a constraint between the reference group
and one or more groups. The setup is identical to the option umbrella,
except for the fact that a rigid constraint is applied instead of a
harmonic potential. Note that this type is not supported in combination
with multiple time stepping.


constant-force

Center of mass pulling using a linear potential and therefore a constant
force. For this option there is no reference position and therefore the
parameters `pull-coord1-init` and `pull-coord1-rate` are not used.


flat-bottom

At distances above `pull-coord1-init` a harmonic potential is applied,
otherwise no potential is applied.


flat-bottom-high

At distances below `pull-coord1-init` a harmonic potential is applied,
otherwise no potential is applied.


external-potential

An external potential that needs to be provided by another module.


pull-coord1-potential-provider

The name of the external module that provides the potential for the case
where `pull-coord1-type=external-potential`.


pull-coord1-geometry

distance

Pull along the vector connecting the two groups. Components can be
selected with `pull-coord1-dim`.


direction

Pull in the direction of `pull-coord1-vec`.


direction-periodic

As `pull-coord1-geometry=direction`, but does not apply periodic box
vector corrections to keep the distance within half the box length. This
is (only) useful for pushing groups apart by more than half the box
length by continuously changing the reference location using a pull
rate. With this geometry the box should not be dynamic (*e.g.* no
pressure scaling) in the pull dimensions and the pull force is not added
to the virial.


direction-relative

As `pull-coord1-geometry=direction`, but the pull vector is the vector
that points from the COM of a third to the COM of a fourth pull group.
This means that 4 groups need to be supplied in `pull-coord1-groups`.
Note that the pull force will give rise to a torque on the pull vector,
which is turn leads to forces perpendicular to the pull vector on the
two groups defining the vector. If you want a pull group to move between
the two groups defining the vector, simply use the union of these two
groups as the reference group.


cylinder

Designed for pulling with respect to a layer where the reference COM is
given by a local cylindrical part of the reference group. The pulling is
in the direction of `pull-coord1-vec`. From the first of the two groups
in `pull-coord1-groups` a cylinder is selected around the axis going
through the COM of the second group with direction `pull-coord1-vec`
with radius `pull-cylinder-r`. Weights of the atoms decrease
continuously to zero as the radial distance goes from 0 to
`pull-cylinder-r` (mass weighting is also used). The radial dependence
gives rise to radial forces on both pull groups. Note that the radius
should be smaller than half the box size. For tilted cylinders they
should be even smaller than half the box size since the distance of an
atom in the reference group from the COM of the pull group has both a
radial and an axial component. This geometry is not supported with
constraint pulling.


angle

Pull along an angle defined by four groups. The angle is defined as the
angle between two vectors: the vector connecting the COM of the first
group to the COM of the second group and the vector connecting the COM
of the third group to the COM of the fourth group.


angle-axis

As `pull-coord1-geometry=angle` but the second vector is given by
`pull-coord1-vec`. Thus, only the two groups that define the first
vector need to be given.


dihedral

Pull along a dihedral angle defined by six groups. These pairwise define
three vectors: the vector connecting the COM of group 1 to the COM of
group 2, the COM of group 3 to the COM of group 4, and the COM of group
5 to the COM group 6. The dihedral angle is then defined as the angle
between two planes: the plane spanned by the the two first vectors and
the plane spanned the two last vectors.


transformation

Transforms other pull coordinates using a mathematical expression
defined by `pull-coord1-expression`. Pull coordinates of lower indices,
and time, can be used as variables to this pull coordinate. Thus, pull
transformation coordinates should have a higher pull coordinate index
than all pull coordinates they transform.


pull-coord1-expression

Mathematical expression to transform pull coordinates of lower indices
to a new one. The pull coordinates are referred to as variables in the
equation so that pull-coord1's value becomes 'x1', pull-coord2 value
becomes 'x2' etc. Time can also be used a variable, becoming 't'. Note
that angular coordinates use units of radians in the expression. The
mathematical expression are evaluated using muParser. Only relevant if
`pull-coord1-geometry=transformation`.


pull-coord1-dx

(1e-9) Size of finite difference to use in numerical derivation of the
pull coordinate with respect to other pull coordinates. The current
implementation uses a simple first order finite difference method to
perform derivation so that f'(x) = (f(x+dx)-f(x))/dx Only relevant if
`pull-coord1-geometry=transformation`.


pull-coord1-groups

The group indices on which this pull coordinate will operate. The number
of group indices required is geometry dependent. The first index is the
reference group and can be 0, in which case an absolute reference of
`pull-coord1-origin` is used. With an absolute reference the system is
no longer translation invariant and one should think about what to do
with the center of mass motion.


pull-coord1-dim

(Y Y Y) Selects the dimensions that this pull coordinate acts on and
that are printed to the output files when `pull-print-components` =
`pull-coord1-start=yes`. With `pull-coord1-geometry` =
`pull-coord1-geometry=distance`, only Cartesian components set to Y
contribute to the distance. Thus setting this to Y Y N results in a
distance in the x/y plane. With other geometries all dimensions with
non-zero entries in `pull-coord1-vec` should be set to Y, the values for
other dimensions only affect the output.


pull-coord1-origin

(0.0 0.0 0.0) The pull reference position for use with an absolute
reference.


pull-coord1-vec

(0.0 0.0 0.0) The pull direction. `gmx grompp` normalizes the vector.


pull-coord1-start

no

do not modify `pull-coord1-init`


yes

add the COM distance of the starting conformation to `pull-coord1-init`


pull-coord1-init

(0.0) \[nm\] or \[deg\] The reference distance or reference angle at
t=0.


pull-coord1-rate

\(0\) \[nm/ps\] or \[deg/ps\] The rate of change of the reference
position or reference angle.


pull-coord1-k

\(0\) \[kJ mol^{-1} nm^{-2}\] or \[kJ mol^{-1}
nm^{-1}\] or \[kJ mol^{-1} rad^{-2}\] or \[kJ
mol^{-1} rad^{-1}\] The force constant. For umbrella
pulling this is the harmonic force constant in kJ mol^{-1}
nm^{-2} (or kJ mol^{-1} rad^{-2} for angles).
For constant force pulling this is the force constant of the linear
potential, and thus the negative (!) of the constant force in kJ
mol^{-1} nm^{-1} (or kJ mol^{-1}
rad^{-1} for angles). Note that for angles the force constant is
expressed in terms of radians (while `pull-coord1-init` and
`pull-coord1-rate` are expressed in degrees).


pull-coord1-kB

(pull-k1) \[kJ mol^{-1} nm^{-2}\] or \[kJ
mol^{-1} nm^{-1}\] or \[kJ mol^{-1}
rad^{-2}\] or \[kJ mol^{-1} rad^{-1}\] As
`pull-coord1-k`, but for state B. This is only used when `free-energy`
is turned on. The force constant is then (1 -lambda) \*
`pull-coord1-k` + lambda \* `pull-coord1-kB`.


### AWH adaptive biasing

awh

no

No biasing.


yes

Adaptively bias a reaction coordinate using the AWH method and estimate
the corresponding PMF. This requires a constant ensemble temperature to
be available. The PMF and other AWH data are written to the energy file
(`edr`) at an interval set by `awh-nstout` and can be extracted with the
`gmx awh` tool. The AWH coordinate can be multidimensional and is
defined by mapping each dimension to a pull coordinate index (and/or up
to one alchemical free lambda state dimension `free-energy`). This is
only allowed if `pull-coord1-type=external-potential` and
`pull-coord1-potential-provider` = `awh` for the concerned pull
coordinate indices. Pull geometry 'direction-periodic' and
transformation coordinates that depend on time are not supported by AWH.


awh-potential

convolved

The applied biasing potential is the convolution of the bias function
and a set of harmonic umbrella potentials (see `awh-potential=umbrella`
below). This results in a smooth potential function and force. The
resolution of the potential is set by the force constant of each
umbrella, see `awh1-dim1-force-constant`. This option is not compatible
with using the free energy lambda state as an AWH reaction coordinate
dimension.


umbrella

The potential bias is applied by controlling the position of a harmonic
potential using Monte-Carlo sampling. The force constant is set with
`awh1-dim1-force-constant`. The umbrella location is sampled using
Monte-Carlo every `awh-nstsample` steps. This is option is required when
using the free energy lambda state as an AWH reaction coordinate
dimension. Apart from that, this option is mainly for comparison and
testing purposes as there are no advantages to using an umbrella.


awh-share-multisim

no

AWH will not share biases across simulations started with `gmx mdrun`
option `-multidir`. The biases will be independent.


yes

With `gmx mdrun` and option `-multidir` the bias and PMF estimates for
biases with `awh1-share-group` \>0 will be shared across simulations
with the biases with the same `awh1-share-group` value. The simulations
should have the same AWH settings for sharing to make sense. `gmx mdrun`
will check whether the simulations are technically compatible for
sharing, but the user should check that bias sharing physically makes
sense.


awh-seed

(-1) Random seed for Monte-Carlo sampling the umbrella position, where
-1 indicates to generate a seed. Only used with
`awh-potential=umbrella`.


awh-nstout

\(100000\) Number of steps between printing AWH data to the energy file,
should be a multiple of `nstenergy`.


awh-nstsample

\(10\) Number of steps between sampling of the coordinate value. This
sampling is the basis for updating the bias and estimating the PMF and
other AWH observables.


awh-nsamples-update

\(100\) The number of coordinate samples used for each AWH update. The
update interval in steps is `awh-nstsample` times this value.


awh-nbias

\(1\) The number of biases, each acting on its own coordinate. The
following options should be specified for each bias although below only
the options for bias number 1 is shown. Options for other bias indices
are obtained by replacing '1' by the bias index.


awh1-error-init

(10.0) \[kJ mol^{-1}\] Estimated initial average error of the
PMF for this bias. This value together with an estimate of the crossing
time, based on the length of the sampling interval and the given
diffusion constant(s) `awh1-dim1-diffusion`, determine the initial
biasing rate. With multiple dimensions, the longest crossing time is
used. The error is obviously not known *a priori*. Only a rough estimate
of `awh1-error-init` is needed however. As a general guideline, leave
`awh1-error-init` to its default value when starting a new simulation.
On the other hand, when there is *a priori* knowledge of the PMF (e.g.
when an initial PMF estimate is provided, see the `awh1-user-data`
option) then `awh1-error-init` should reflect that knowledge.


awh1-growth

Each bias keeps a reference weight histogram for the coordinate samples.
Its size sets the magnitude of the bias function and free energy
estimate updates (few samples corresponds to large updates and vice
versa). Thus, its growth rate sets the maximum convergence rate.

exp-linear

By default, there is an initial stage in which the histogram grows close
to exponentially (but slower than the sampling rate). In the final stage
that follows, the growth rate is linear and equal to the sampling rate
(set by `awh-nstsample`). The initial stage is typically necessary for
efficient convergence when starting a new simulation where high free
energy barriers have not yet been flattened by the bias.


linear

As `awh1-growth=exp-linear` but skip the initial stage. This may be
useful if there is *a priori* knowledge (see `awh1-error-init`) which
eliminates the need for an initial stage. This is also the setting
compatible with `awh1-target=local-boltzmann`.


awh1-growth-factor

\(2\) \[\] The growth factor $\gamma$ during the exponential phase
with `awh1-growth=exp-linear`. Should be larger than 1.


awh1-equilibrate-histogram

yes

Before entering the initial stage (see `awh1-growth=exp-linear`), make
sure the histogram of sampled weights is following the target
distribution closely enough (specifically, at least 80% of the target
region needs to have a local relative error of less than
`awh1-histogram-tolerance`). This option is particularly important when
`awh1-share-group` \> 0 and the initial configurations poorly represent
the target distribution.


no

Do not equilibrate the histogram.


awh1-histogram-tolerance

(0.3) \[\] The relative tolerance for the histogram of sampled weigths,
used with `awh1-equilibrate-histogram=yes`. The value of 0.3 is
sufficiently large to not slow down the convergence, even when a single
walker is used.


awh1-target

constant

The bias is tuned towards a constant (uniform) coordinate distribution
in the defined sampling interval (defined by \[`awh1-dim1-start`,
`awh1-dim1-end`\]).


cutoff

Similar to `awh1-target=constant`, but the target distribution is
proportional to 1/(1 + exp(F - `awh1-target=cutoff`)), where F is the
free energy relative to the estimated global minimum. This provides a
smooth switch of a flat target distribution in regions with free energy
lower than the cut-off to a Boltzmann distribution in regions with free
energy higher than the cut-off.


boltzmann

The target distribution is a Boltzmann distribtution with a scaled beta
(inverse temperature) factor given by `awh1-target-beta-scaling`.
*E.g.*, a value of 0.1 would give the same coordinate distribution as
sampling with a simulation temperature scaled by 10.


local-boltzmann

Same target distribution and use of `awh1-target-beta-scaling` but the
convergence towards the target distribution is inherently local *i.e.*,
the rate of change of the bias only depends on the local sampling. This
local convergence property is only compatible with `awh1-growth=linear`,
since for `awh1-growth=exp-linear` histograms are globally rescaled in
the initial stage.


awh1-target-beta-scaling

\(0\) For `awh1-target=boltzmann` and `awh1-target=local-boltzmann` it
is the unitless beta scaling factor taking values in (0,1).


awh1-target-cutoff

\(0\) \[kJ mol^{-1}\] For `awh1-target=cutoff` this is the
cutoff, should be \> 0.


awh1-user-data

no

Initialize the PMF and target distribution with default values.


yes

Initialize the PMF and target distribution with user provided data. For
`awh-nbias` = 1, `gmx mdrun` will expect a file `awhinit.xvg` to be
present in the run directory. For multiple biases, `gmx mdrun` expects
files `awhinit1.xvg`, `awhinit2.xvg`, etc. The file name can be changed
with the `-awh` option. The first `awh1-ndim` columns of each input file
should contain the coordinate values, such that each row defines a point
in coordinate space. Column `awh1-ndim` + 1 should contain the PMF value
(in kT) for each point. The target distribution column can either follow
the PMF (column `awh1-ndim` + 2) or be in the same column as written by
`gmx awh`.


awh1-share-group

0

Do not share the bias.


positive

Share the bias and PMF estimates between simulations. This currently
only works between biases with the same index. Note that currently
sharing within a single simulation is not supported. The bias will be
shared across simulations that specify the same value for
`awh1-share-group`. To enable this, use `awh-share-multisim=yes` and the
`gmx mdrun` option `-multidir`. Sharing may increase convergence
initially, although the starting configurations can be critical,
especially when sharing between many biases. N.b., multiple walkers
sharing a degenerate reaction coordinate may have problems overlapping
their sampling, possibly making it difficult to cover the sampling
interval.


awh1-target-metric-scaling

no

Do not scale the target distribution based on the AWH friction metric.


yes

Scale the target distribution based on the AWH friction metric. Regions
with high friction (long autocorrelation times) will be sampled more.
The diffusion metric is the inverse of the friction metric. This scaling
can be used with any `awh1-target` type and is applied after
user-provided target distribution modifications (`awh1-user-data`), if
any. If `awh1-growth=exp-linear`, the target distribution scaling starts
after leaving the initial phase.


awh1-target-metric-scaling-limit

\(10\) The upper limit of scaling, relative to the average, when
`awh1-target-metric-scaling` is enabled. The lower limit will be the
inverse of this value. This upper limit should be \> 1.


awh1-ndim

\(1\) \[integer\] Number of dimensions of the coordinate, each dimension
maps to 1 pull coordinate. The following options should be specified for
each such dimension. Below only the options for dimension number 1 is
shown. Options for other dimension indices are obtained by replacing '1'
by the dimension index.


awh1-dim1-coord-provider

pull

The pull module is providing the reaction coordinate for this dimension.
With multiple time-stepping, AWH and pull should be in the same MTS
level.


fep-lambda

The free energy `free-energy` lambda state is the reaction coordinate
for this dimension. The lambda states to use are specified by
`fep-lambdas`, `vdw-lambdas`, `coul-lambdas` etc. This is not compatible
with delta-lambda. It also requires `calc-lambda-neighbors=-1`. With
multiple time-stepping, AWH should be in the slow level. This option
requires `awh-potential=umbrella`.


awh1-dim1-coord-index

\(1\) Index of the pull coordinate defining this coordinate dimension.


awh1-dim1-force-constant

\(0\) \[kJ mol^{-1} nm^{-2}\] or \[kJ mol^{-1}
rad^{-2}\] Force constant for the (convolved) umbrella
potential(s) along this coordinate dimension.


awh1-dim1-start

(0.0) \[nm\] or \[deg\] Start value of the sampling interval along this
dimension. The range of allowed values depends on the relevant pull
geometry (see `pull-coord1-geometry`). For dihedral geometries
`awh1-dim1-start` greater than `awh1-dim1-end` is allowed. The interval
will then wrap around from +period/2 to -period/2. For the direction
geometry, the dimension is made periodic when the direction is along a
box vector and covers more than 95% of the box length. Note that one
should not apply pressure coupling along a periodic dimension.


awh1-dim1-end

(0.0) \[nm\] or \[deg\] End value defining the sampling interval
together with `awh1-dim1-start`.


awh1-dim1-diffusion

(10^{-5}) \[nm^{2}/ps\], \[rad^{2}/ps\] or
\[ps^{-1}\] Estimated diffusion constant for this coordinate
dimension determining the initial biasing rate. This needs only be a
rough estimate and should not critically affect the results unless it is
set to something very low, leading to slow convergence, or very high,
forcing the system far from equilibrium. Not setting this value
explicitly generates a warning.


awh1-dim1-cover-diameter

(0.0) \[nm\] or \[deg\] Diameter that needs to be sampled by a single
simulation around a coordinate value before the point is considered
covered in the initial stage (see `awh1-growth=exp-linear`). A value \>
0 ensures that for each covering there is a continuous transition of
this diameter across each coordinate value. This is trivially true for
independent simulations but not for for multiple bias-sharing
simulations (`awh1-share-group`\>0). For a diameter = 0, covering occurs
as soon as the simulations have sampled the whole interval, which for
many sharing simulations does not guarantee transitions across free
energy barriers. On the other hand, when the diameter \>= the sampling
interval length, covering occurs when a single simulation has
independently sampled the whole interval.


### Enforced rotation

These `mdp` parameters can be used enforce the rotation of a group of
atoms, e.g. a protein subunit. The [reference manual]() describes in
detail 13 different potentials that can be used to achieve such a
rotation.

rotation

no

No enforced rotation will be applied. All enforced rotation options will
be ignored (and if present in the `mdp` file, they unfortunately
generate warnings).


yes

Apply the rotation potential specified by `rot-type0` to the group of
atoms given under the `rot-group0` option.


rot-ngroups

\(1\) Number of rotation groups.


rot-group0

Name of rotation group 0 in the index file.


rot-type0

(iso) Type of rotation potential that is applied to rotation group 0.
Can be of of the following: `iso`, `iso-pf`, `pm`, `pm-pf`, `rm`,
`rm-pf`, `rm2`, `rm2-pf`, `flex`, `flex-t`, `flex2`, or `flex2-t`.


rot-massw0

(no) Use mass weighted rotation group positions.


rot-vec0

(1.0 0.0 0.0) Rotation vector, will get normalized.


rot-pivot0

(0.0 0.0 0.0) \[nm\] Pivot point for the potentials `iso`, `pm`, `rm`,
and `rm2`.


rot-rate0

\(0\) \[degree ps^{-1}\] Reference rotation rate of group 0.


rot-k0

\(0\) \[kJ mol^{-1} nm^{-2}\] Force constant for group
0.


rot-slab-dist0

(1.5) \[nm\] Slab distance, if a flexible axis rotation type was chosen.


rot-min-gauss0

(0.001) Minimum value (cutoff) of Gaussian function for the force to be
evaluated (for the flexible axis potentials).


rot-eps0

(0.0001) \[nm^{2}\] Value of additive constant epsilon for
`rm2*` and `flex2*` potentials.


rot-fit-method0

(rmsd) Fitting method when determining the actual angle of a rotation
group (can be one of `rmsd`, `norm`, or `potential`).


rot-potfit-nsteps0

\(21\) For fit type `potential`, the number of angular positions around
the reference angle for which the rotation potential is evaluated.


rot-potfit-step0

(0.25) For fit type `potential`, the distance in degrees between two
angular positions.


rot-nstrout

\(100\) Output interval (in steps) for the angle of the rotation group,
as well as for the torque and the rotation potential energy.


rot-nstsout

\(1000\) Output interval (in steps) for per-slab data of the flexible
axis potentials, i.e. angles, torques and slab centers.


### NMR refinement

disre

no

ignore distance restraint information in topology file


simple

simple (per-molecule) distance restraints.


ensemble

distance restraints over an ensemble of molecules in one simulation box.
Normally, one would perform ensemble averaging over multiple
simulations, using `mdrun -multidir`. The environment variable
`GMX_DISRE_ENSEMBLE_SIZE` sets the number of systems within each
ensemble (usually equal to the number of directories supplied to
`mdrun -multidir`).


disre-weighting

equal

divide the restraint force equally over all atom pairs in the restraint


conservative

the forces are the derivative of the restraint potential, this results
in an weighting of the atom pairs to the reciprocal seventh power of the
displacement. The forces are conservative when `disre-tau` is zero.


disre-mixed

no

the violation used in the calculation of the restraint force is the
time-averaged violation


yes

the violation used in the calculation of the restraint force is the
square root of the product of the time-averaged violation and the
instantaneous violation


disre-fc

\(1000\) \[kJ mol^{-1} nm^{-2}\] force constant for
distance restraints, which is multiplied by a (possibly) different
factor for each restraint given in the `fac` column of the interaction
in the topology file.


disre-tau

\(0\) \[ps\] time constant for distance restraints running average. A
value of zero turns off time averaging.


nstdisreout

\(100\) \[steps\] period between steps when the running time-averaged
and instantaneous distances of all atom pairs involved in restraints are
written to the energy file (can make the energy file very large)


orire

no

ignore orientation restraint information in topology file


yes

use orientation restraints, ensemble averaging can be performed with
`mdrun -multidir`


orire-fc

\(0\) \[kJ mol^{-1}\] force constant for orientation restraints,
which is multiplied by a (possibly) different weight factor for each
restraint, can be set to zero to obtain the orientations from a free
simulation


orire-tau

\(0\) \[ps\] time constant for orientation restraints running average. A
value of zero turns off time averaging.


orire-fitgrp

fit group for orientation restraining. This group of atoms is used to
determine the rotation **R** of the system with respect to the reference
orientation. The reference orientation is the starting conformation of
the first subsystem. For a protein, backbone is a reasonable choice


nstorireout

\(100\) \[steps\] period between steps when the running time-averaged
and instantaneous orientations for all restraints, and the molecular
order tensor are written to the energy file (can make the energy file
very large)


### Free energy calculations

free-energy

no

Only use topology A.


yes

Interpolate between topology A (lambda=0) to topology B (lambda=1) and
write the derivative of the Hamiltonian with respect to lambda (as
specified with `dhdl-derivatives`), or the Hamiltonian differences with
respect to other "foreign" lambda values (as specified with
`calc-lambda-neighbors`) to the energy file and/or to `dhdl.xvg`, where
they can be processed by, for example `gmx bar`. The potentials,
bond-lengths and angles are interpolated linearly as described in the
manual. When `sc-alpha` is larger than zero, soft-core potentials are
used for the LJ and Coulomb interactions.


expanded

Turns on expanded ensemble simulation, where the alchemical state
becomes a dynamic variable, allowing jumping between different
Hamiltonians. See the expanded ensemble options for controlling how
expanded ensemble simulations are performed. The different Hamiltonians
used in expanded ensemble simulations are defined by the other free
energy options.


init-lambda

(-1) starting value for lambda (float). Generally, this should only be
used with slow growth (*i.e.* nonzero `delta-lambda`). In other cases,
`init-lambda-state` should be specified instead. If a lambda vector is
given, `init-lambda` is used to interpolate the vector instead of
setting lambda directly. Must be greater than or equal to 0.


delta-lambda

\(0\) increment per time step for lambda


init-lambda-state

(-1) starting value for the lambda state (integer). Specifies which
column of the lambda vector (`coul-lambdas`, `vdw-lambdas`,
`bonded-lambdas`, `restraint-lambdas`, `mass-lambdas`,
`temperature-lambdas`, `fep-lambdas`) should be used. This is a
zero-based index: `init-lambda-state=0` means the first column, and so
on.


fep-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. Free energy differences between different
lambda values can then be determined with `gmx bar`. `fep-lambdas` is
different from the other -lambdas keywords because all components of the
lambda vector that are not specified will use `fep-lambdas`.


coul-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. If soft-core potentials are used, values must
be between 0 and 1. Only the electrostatic interactions are controlled
with this component of the lambda vector (and only if the lambda=0 and
lambda=1 states have differing electrostatic interactions).


vdw-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. If soft-core potentials are used, values must
be between 0 and 1. Only the van der Waals interactions are controlled
with this component of the lambda vector.


bonded-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. Only the bonded interactions are controlled
with this component of the lambda vector.


restraint-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. Only the restraint interactions: dihedral
restraints, and the pull code restraints are controlled with this
component of the lambda vector.


mass-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. Only the particle masses are controlled with
this component of the lambda vector.


temperature-lambdas

\[array\] Zero, one or more lambda values for which Delta H values will
be determined and written to dhdl.xvg every `nstdhdl` steps. Values must
be greater than or equal to 0; values greater than 1 are allowed but
should be used carefully. Only the temperatures are controlled with this
component of the lambda vector. Note that these lambdas should not be
used for replica exchange, only for simulated tempering.


calc-lambda-neighbors

\(1\) Controls the number of lambda values for which Delta H values will
be calculated and written out, if `init-lambda-state` has been set.
These lambda values are referred to as "foreign" lambdas. A positive
value will limit the number of lambda points calculated to only the nth
neighbors of `init-lambda-state`: for example, if `init-lambda-state` is
5 and this parameter has a value of 2, energies for lambda points 3-7
will be calculated and writen out. A value of -1 means all lambda points
will be written out. For normal BAR such as with `gmx bar`, a value of 1
is sufficient, while for MBAR -1 should be used.


sc-function

(beutler)

beutler


Beutler *et al.* soft-core function

gapsys


Gapsys *et al.* soft-core function


sc-alpha

\(0\) for `sc-function=beutler` the
soft-core alpha parameter, a value of 0 results in linear interpolation
of the LJ and Coulomb interactions. Used only with
`sc-function=beutler`


sc-r-power

\(6\) power 6 for the radial term in the soft-core equation. Used only
with `sc-function=beutler`


sc-coul

(no) Whether to apply the soft-core free energy interaction
transformation to the Coulombic interaction of a molecule. Default is
no, as it is generally more efficient to turn off the Coulombic
interactions linearly before turning off the van der Waals interactions.
Note that it is only taken into account when there are multiple lambda
components, and you can still turn off soft-core interactions by setting
`sc-alpha` to 0. Used only with
`sc-function=beutler`


sc-power

\(1\) the power for lambda in the soft-core function, only the values 1
and 2 are supported. Used only with
`sc-function=beutler`


sc-sigma

(0.3) \[nm\] for `sc-function=beutler` the
soft-core sigma for particles which have a C6 or C12 parameter equal to
zero or a sigma smaller than `sc-sigma`. Used only with
`sc-function=beutler`


sc-gapsys-scale-linpoint-lj

(0.85) for `sc-function=gapsys` it is the
unitless alphaLJ parameter. It controls the softness of the van der
Waals interactions by scaling the point for linearizing the vdw force.
Setting it to 0 will result in the standard hard-core van der Waals
interactions. Used only with
`sc-function=gapsys`


sc-gapsys-scale-linpoint-q

(0.3) \[nm/e^2\] For `sc-function=gapsys`
the alphaQ parameter with a default value of 0.3. It controls the
softness of the Coulombic interactions. Setting it to 0 will result in
the standard hard-core Coulombic interactions. Used only with
`sc-function=gapsys`


sc-gapsys-sigma-lj

(0.3) \[nm\] for `sc-function=gapsys` the
soft-core sigma for particles which have a C6 or C12 parameter equal to
zero. Used only with `sc-function=gapsys`


couple-moltype

Here one can supply a molecule type (as defined in the topology) for
calculating solvation or coupling free energies. There is a special
option `system` that couples all molecule types in the system. This can
be useful for equilibrating a system starting from (nearly) random
coordinates. `free-energy` has to be turned on. The Van der Waals
interactions and/or charges in this molecule type can be turned on or
off between lambda=0 and lambda=1, depending on the settings of
`couple-lambda0` and `couple-lambda1`. If you want to decouple one of
several copies of a molecule, you need to copy and rename the molecule
definition in the topology.


couple-lambda0

vdw-q

all interactions are on at lambda=0


vdw

the charges are zero (no Coulomb interactions) at lambda=0


q

the Van der Waals interactions are turned off at lambda=0; soft-core
interactions will be required to avoid singularities


none

the Van der Waals interactions are turned off and the charges are zero
at lambda=0; soft-core interactions will be required to avoid
singularities.


couple-lambda1

analogous to `couple-lambda0`, but for lambda=1


couple-intramol

no

All intra-molecular non-bonded interactions for moleculetype
`couple-moltype` are replaced by exclusions and explicit pair
interactions. In this manner the decoupled state of the molecule
corresponds to the proper vacuum state without periodicity effects.


yes

The intra-molecular Van der Waals and Coulomb interactions are also
turned on/off. This can be useful for partitioning free-energies of
relatively large molecules, where the intra-molecular non-bonded
interactions might lead to kinetically trapped vacuum conformations. The
1-4 pair interactions are not turned off.


nstdhdl

\(100\) the interval for writing dH/dlambda and possibly Delta H to
dhdl.xvg, 0 means no ouput, should be a multiple of `nstcalcenergy`.


dhdl-derivatives

(yes)

If yes (the default), the derivatives of the Hamiltonian with respect to
lambda at each `nstdhdl` step are written out. These values are needed
for interpolation of linear energy differences with `gmx bar` (although
the same can also be achieved with the right `calc-lambda-neighbors`
setting, that may not be as flexible), or with thermodynamic integration


dhdl-print-energy

(no)

Include either the total or the potential energy in the dhdl file.
Options are 'no', 'potential', or 'total'. This information is needed
for later free energy analysis if the states of interest are at
different temperatures. If all states are at the same temperature, this
information is not needed. 'potential' is useful in case one is using
`mdrun -rerun` to generate the `dhdl.xvg` file. When rerunning from an
existing trajectory, the kinetic energy will often not be correct, and
thus one must compute the residual free energy from the potential alone,
with the kinetic energy component computed analytically.


separate-dhdl-file

yes

The free energy values that are calculated (as specified with
`calc-lambda-neighbors` and `dhdl-derivatives` settings) are written out
to a separate file, with the default name `dhdl.xvg`. This file can be
used directly with `gmx
bar`.


no

The free energy values are written out to the energy output file
(`ener.edr`, in accumulated blocks at every `nstenergy` steps), where
they can be extracted with `gmx energy` or used directly with `gmx bar`.


dh-hist-size

\(0\) If nonzero, specifies the size of the histogram into which the
Delta H values (specified with `calc-lambda-neighbors`) and the
derivative dH/dl values are binned, and written to ener.edr. This can be
used to save disk space while calculating free energy differences. One
histogram gets written for each foreign lambda and two for the dH/dl, at
every `nstenergy` step. Be aware that incorrect histogram settings (too
small size or too wide bins) can introduce errors. Do not use histograms
unless you are certain you need it.


dh-hist-spacing

(0.1) Specifies the bin width of the histograms, in energy units. Used
in conjunction with `dh-hist-size`. This size limits the accuracy with
which free energies can be calculated. Do not use histograms unless you
are certain you need it.


### Expanded Ensemble calculations

nstexpanded

The number of integration steps beween attempted moves changing the
system Hamiltonian in expanded ensemble simulations. Must be a multiple
of `nstcalcenergy`, but can be greater or less than `nstdhdl`.


lmc-stats

no

No Monte Carlo in state space is performed.


metropolis-transition

Uses the Metropolis weights to update the expanded ensemble weight of
each state. Min{1,exp(-(beta_new u_new - beta_old u_old)}


barker-transition

Uses the Barker transition critera to update the expanded ensemble
weight of each state i, defined by exp(-beta_new u_new)/(exp(-beta_new
u_new)+exp(-beta_old u_old))


wang-landau

Uses the Wang-Landau algorithm (in state space, not energy space) to
update the expanded ensemble weights.


min-variance

Uses the minimum variance updating method of Escobedo et al. to update
the expanded ensemble weights. Weights will not be the free energies,
but will rather emphasize states that need more sampling to give even
uncertainty.


lmc-mc-move

no

No Monte Carlo in state space is performed.


metropolis-transition

Randomly chooses a new state up or down, then uses the Metropolis
criteria to decide whether to accept or reject: Min{1,exp(-(beta_new
u_new - beta_old u_old)}


barker-transition

Randomly chooses a new state up or down, then uses the Barker transition
criteria to decide whether to accept or reject: exp(-beta_new
u_new)/(exp(-beta_new u_new)+exp(-beta_old u_old))


gibbs

Uses the conditional weights of the state given the coordinate
(exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move
to.


metropolized-gibbs

Uses the conditional weights of the state given the coordinate
(exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move
to, EXCLUDING the current state, then uses a rejection step to ensure
detailed balance. Always more efficient that Gibbs, though only
marginally so in many situations, such as when only the nearest
neighbors have decent phase space overlap.


lmc-seed

(-1) random seed to use for Monte Carlo moves in state space. When
`lmc-seed` is set to -1, a pseudo random seed is us


mc-temperature

Temperature used for acceptance/rejection for Monte Carlo moves. If not
specified, the temperature of the simulation specified in the first
group of `ref-t` is used.


wl-ratio

(0.8) The cutoff for the histogram of state occupancies to be reset, and
the free energy incrementor to be changed from delta to delta \*
`wl-scale`. If we define the Nratio = (number of samples at each
histogram) / (average number of samples at each histogram). `wl-ratio`
of 0.8 means that means that the histogram is only considered flat if
all Nratio \> 0.8 AND simultaneously all 1/Nratio \> 0.8.


wl-scale

(0.8) Each time the histogram is considered flat, then the current value
of the Wang-Landau incrementor for the free energies is multiplied by
`wl-scale`. Value must be between 0 and 1.


init-wl-delta

(1.0) The initial value of the Wang-Landau incrementor in kT. Some value
near 1 kT is usually most efficient, though sometimes a value of 2-3 in
units of kT works better if the free energy differences are large.


wl-oneovert

(no) Set Wang-Landau incrementor to scale with 1/(simulation time) in
the large sample limit. There is significant evidence that the standard
Wang-Landau algorithms in state space presented here result in free
energies getting 'burned in' to incorrect values that depend on the
initial state. when `wl-oneovert` is true, then when the incrementor
becomes less than 1/N, where N is the number of samples collected (and
thus proportional to the data collection time, hence '1 over t'), then
the Wang-Lambda incrementor is set to 1/N, decreasing every step. Once
this occurs, `wl-ratio` is ignored, but the weights will still stop
updating when the equilibration criteria set in `lmc-weights-equil` is
achieved.


lmc-repeats

\(1\) Controls the number of times that each Monte Carlo swap type is
performed each iteration. In the limit of large numbers of Monte Carlo
repeats, then all methods converge to Gibbs sampling. The value will
generally not need to be different from 1.


lmc-gibbsdelta

(-1) Limit Gibbs sampling to selected numbers of neighboring states. For
Gibbs sampling, it is sometimes inefficient to perform Gibbs sampling
over all of the states that are defined. A positive value of
`lmc-gibbsdelta` means that only states plus or minus `lmc-gibbsdelta`
are considered in exchanges up and down. A value of -1 means that all
states are considered. For less than 100 states, it is probably not that
expensive to include all states.


lmc-forced-nstart

\(0\) Force initial state space sampling to generate weights. In order
to come up with reasonable initial weights, this setting allows the
simulation to drive from the initial to the final lambda state, with
`lmc-forced-nstart` steps at each state before moving on to the next
lambda state. If `lmc-forced-nstart` is sufficiently long (thousands of
steps, perhaps), then the weights will be close to correct. However, in
most cases, it is probably better to simply run the standard weight
equilibration algorithms.


nst-transition-matrix

(-1) Interval of outputting the expanded ensemble transition matrix. A
negative number means it will only be printed at the end of the
simulation.


symmetrized-transition-matrix

(no) Whether to symmetrize the empirical transition matrix. In the
infinite limit the matrix will be symmetric, but will diverge with
statistical noise for short timescales. Forced symmetrization, by using
the matrix T_sym = 1/2 (T + transpose(T)), removes problems like the
existence of (small magnitude) negative eigenvalues.


mininum-var-min

\(100\) The min-variance strategy (option of `lmc-stats` is only valid
for larger number of samples, and can get stuck if too few samples are
used at each state. `mininum-var-min` is the minimum number of samples
that each state that are allowed before the min-variance strategy is
activated if selected.


init-lambda-weights

The initial weights (free energies) used for the expanded ensemble
states. Default is a vector of zero weights. format is similar to the
lambda vector settings in `fep-lambdas`, except the weights can be any
floating point number. Units are kT. Its length must match the lambda
vector lengths.


init-wl-histogram-counts

The initial counts used for the Wang-Landau histogram of visiting
expanded ensemble states. The flatness of this histogram is used to
decide whether to decrement the histogram-building incrementor. This
option is only generally useful if continuing a shorter simulation from
a previous one, as the smaller the incrementor gets, the longer it takes
for the histogram to become flat, often longer than a short simulation
takes, requiring the histogram population to be carried over from the
previous simulation. The default is a vector of zeros. The format is
similar to the lambda vector settings in `fep-lambdas`. The value can be
a floating point number or an integer, as some methods increment
multiple histogram bins at the same time with fractional weights. Its
length must match the lambda vector lengths.


init-lambda-counts

The initial counts used for the number of times each expanded ensemble
state is visited states. Several algorithms set by `lmc-weights-equil`
use various functions of the number of visits to each state states to
decide whether to switch to different phases of weight determination.
These include `number-all-lambda` which requires the mumber of times
each lambda state is visited to be equal to or greater than this number,
`number-samples`, which requires the total number of visits to all
lambda states to be greater than or equal to this, and `count-ratio`,
which requires the number of states visited at each state to be within a
given ratio of equal visitation. This option is only generally useful if
continuing a shorter simulation from a previous one, as most methods
will reach the triggering conditions with relatively low number of
samples collected. The default is a vector of zeros. The format is
similar to the lambda vector settings in `fep-lambdas`. Unlike
`init-wl-histogram`, the value can only be an integer. Its length must
match the lambda vector lengths.


lmc-weights-equil

no

Expanded ensemble weights continue to be updated throughout the
simulation.


yes

The input expanded ensemble weights are treated as equilibrated, and are
not updated throughout the simulation.


wl-delta

Expanded ensemble weight updating is stopped when the Wang-Landau
incrementor falls below this value.


number-all-lambda

Expanded ensemble weight updating is stopped when the number of samples
at all of the lambda states is greater than this value.


number-steps

Expanded ensemble weight updating is stopped when the number of steps is
greater than the level specified by this value.


number-samples

Expanded ensemble weight updating is stopped when the number of total
samples across all lambda states is greater than the level specified by
this value.


count-ratio

Expanded ensemble weight updating is stopped when the ratio of samples
at the least sampled lambda state and most sampled lambda state greater
than this value.


simulated-tempering

(no) Turn simulated tempering on or off. Simulated tempering is
implemented as expanded ensemble sampling with different temperatures
instead of different Hamiltonians.


sim-temp-low

\(300\) \[K\] Low temperature for simulated tempering.


sim-temp-high

\(300\) \[K\] High temperature for simulated tempering.


simulated-tempering-scaling

Controls the way that the temperatures at intermediate lambdas are
calculated from the `temperature-lambdas` part of the lambda vector.

linear

Linearly interpolates the temperatures using the values of
`temperature-lambdas`, *i.e.* if `sim-temp-low` =300, `sim-temp-high`
=400, then lambda=0.5 correspond to a temperature of 350. A nonlinear
set of temperatures can always be implemented with uneven spacing in
lambda.


geometric

Interpolates temperatures geometrically between `sim-temp-low` and
`sim-temp-high`. The i:th state has temperature `sim-temp-low` \*
(`sim-temp-high` / `sim-temp-low`) raised to the power of
(i/(ntemps-1)). This should give roughly equal exchange for constant
heat capacity, though of course things simulations that involve protein
folding have very high heat capacity peaks.


exponential

Interpolates temperatures exponentially between `sim-temp-low` and
`sim-temp-high`. The i:th state has temperature `sim-temp-low` +
(`sim-temp-high` -`sim-temp-low`)\*((exp(`temperature-lambdas`
(i))-1)/(exp(1.0)-i)).


### Non-equilibrium MD

acc-grps

groups for constant acceleration (*e.g.* `Protein Sol`) all atoms in
groups Protein and Sol will experience constant acceleration as
specified in the `accelerate` line. Note that the kinetic energy of the
center of mass of accelerated groups contributes to the kinetic energy
and temperature of the system. If this is not desired, make each
accelerate group also a separate temperature coupling group.


accelerate

\(0\) \[nm ps^{-2}\] acceleration for `acc-grps`; x, y and z for
each group (*e.g.* `0.1 0.0 0.0 -0.1 0.0 0.0` means that first group has
constant acceleration of 0.1 nm ps^{-2} in X direction, second
group the opposite).


freezegrps

Groups that are to be frozen (*i.e.* their X, Y, and/or Z position will
not be updated; *e.g.* `Lipid SOL`). `freezedim` specifies for which
dimension(s) the freezing applies. Note that the virial and pressure are
usually not meaningful when frozen atoms are present. Note that
coordinates of frozen atoms are not scaled by pressure-coupling
algorithms.


freezedim

dimensions for which groups in `freezegrps` should be frozen, specify
`Y` or `N` for X, Y and Z and for each group (*e.g.* `Y Y N N N N` means
that particles in the first group can move only in Z direction. The
particles in the second group can move in any direction).


cos-acceleration

\(0\) \[nm ps^{-2}\] the amplitude of the acceleration profile
for calculating the viscosity. The acceleration is in the X-direction
and the magnitude is `cos-acceleration` cos(2 pi z/boxheight). Two terms
are added to the energy file: the amplitude of the velocity profile and
1/viscosity.


deform

(0 0 0 0 0 0) \[nm ps^{-1}\] The velocities of deformation for
the box elements: a(x) b(y) c(z) b(x) c(x) c(y). Each step the box
elements for which `deform` is non-zero are calculated as:
box(ts)+(t-ts)\*deform, off-diagonal elements are corrected for
periodicity. The time ts is set to t at the first step and at steps at
which x and v are written to trajectory to ensure exact restarts.
Deformation can be used together with semiisotropic or anisotropic
pressure coupling when the appropriate compressibilities are set to
zero. The diagonal elements can be used to strain a solid. The
off-diagonal elements can be used to shear a solid or a liquid. Note
that the atom positions are not affected directly by this option.
Instead, the deform option only modifies the velocities of particles
that are shifted by a periodic box vector such that their new velocities
match the virtual velocity flow field corresponding to the box
deformation. As the deform option never accelerates the remaining
particles in the system, the matching velocity flow field should be set
up at the beginning of the simulation to make the particles follow the
deformation. This can be done with the `deform-init-flow` option. The
flow field is removed from the kinetic energy by `gmx mdrun` so the
actual temperature and pressure of the system are reported.


deform-init-flow

no

Do not modify the velocities. Only use this option when the velocities
of the atoms in the initial configuration already obey the flow field.


yes

When the `deform` option is active, add a velocity profile corresponding
to the box deformation to the initial velocities. This is done after
computing observables from the initial state such as the initial
temperature.


### Electric fields

electric-field-x


electric-field-y


electric-field-z

Here you can specify an electric field that optionally can be
alternating and pulsed. The general expression for the field has the
form of a gaussian laser pulse:

$$
E(t) = E_0 \exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right]\cos\left[\omega (t-t_0)\right]
$$

For example, the four parameters for direction x are set in the fields
of `electric-field-x` (and similar for `electric-field-y` and
`electric-field-z`) like

`electric-field-x  = E0 omega t0 sigma`

with units (respectively) V nm^{-1}, ps^{-1}, ps, ps.

In the special case that `sigma = 0`, the exponential term is omitted
and only the cosine term is used. In this case, `t0` must be set to 0.
If also `omega = 0` a static electric field is applied.

Read more at `electric fields` and in ref. `146 <refCaleman2008a>`.


### Mixed quantum/classical molecular dynamics

QMMM-grps

groups to be described at the QM level for MiMiC QM/MM


QMMM

no

QM/MM is no longer supported via these .mdp options. For MiMic, use no
here.


### Computational Electrophysiology

Use these options to switch on and control ion/water position exchanges
in "Computational Electrophysiology" simulation setups. (See the
[reference manual]() for details).

swapcoords

no

Do not enable ion/water position exchanges.


X ; Y ; Z

Allow for ion/water position exchanges along the chosen direction. In a
typical setup with the membranes parallel to the x-y plane, ion/water
pairs need to be exchanged in Z direction to sustain the requested ion
concentrations in the compartments.


swap-frequency

\(1\) The swap attempt frequency, i.e. every how many time steps the ion
counts per compartment are determined and exchanges made if necessary.
Normally, it is not necessary to check at every time step. For typical
Computational Electrophysiology setups, a value of about 100 is
sufficient and yields a negligible performance impact.


split-group0

Name of the index group of the membrane-embedded part of channel \#0.
The center of mass of these atoms defines one of the compartment
boundaries and should be chosen such that it is near the center of the
membrane.


split-group1

Defines the position of the other compartment boundary.


massw-split0

(no) Defines whether or not mass-weighting is used to calculate the
split group center.

no

Use the geometrical center.


yes

Use the center of mass.


massw-split1

(no) As above, but for `split-group1`.


solvent-group

Name of the index group of solvent molecules.


coupl-steps

\(10\) Average the number of ions per compartment over these many swap
attempt steps. This can be used to prevent that ions near a compartment
boundary (diffusing through a channel, e.g.) lead to unwanted back and
forth swaps.


iontypes

\(1\) The number of different ion types to be controlled. These are
during the simulation exchanged with solvent molecules to reach the
desired reference numbers.


iontype0-name

Name of the first ion type.


iontype0-in-A

(-1) Requested (=reference) number of ions of type 0 in compartment A.
The default value of -1 means: use the number of ions as found in time
step 0 as reference value.


iontype0-in-B

(-1) Reference number of ions of type 0 for compartment B.


bulk-offsetA

(0.0) Offset of the first swap layer from the compartment A midplane. By
default (i.e. bulk offset = 0.0), ion/water exchanges happen between
layers at maximum distance (= bulk concentration) to the split group
layers. However, an offset b (-1.0 \< b \< +1.0) can be specified to
offset the bulk layer from the middle at 0.0 towards one of the
compartment-partitioning layers (at +/- 1.0).


bulk-offsetB

(0.0) Offset of the other swap layer from the compartment B midplane.


threshold

\(1\) Only swap ions if threshold difference to requested count is
reached.


cyl0-r

(2.0) \[nm\] Radius of the split cylinder \#0. Two split cylinders
(mimicking the channel pores) can optionally be defined relative to the
center of the split group. With the help of these cylinders it can be
counted which ions have passed which channel. The split cylinder
definition has no impact on whether or not ion/water swaps are done.


cyl0-up

(1.0) \[nm\] Upper extension of the split cylinder \#0.


cyl0-down

(1.0) \[nm\] Lower extension of the split cylinder \#0.


cyl1-r

(2.0) \[nm\] Radius of the split cylinder \#1.


cyl1-up

(1.0) \[nm\] Upper extension of the split cylinder \#1.


cyl1-down

(1.0) \[nm\] Lower extension of the split cylinder \#1.


### Density-guided simulations

These options enable and control the calculation and application of
additional forces that are derived from three-dimensional densities,
e.g., from cryo electron-microscopy experiments. (See the [reference
manual]() for details)

density-guided-simulation-active

(no) Activate density-guided simulations.


density-guided-simulation-group

(protein) The atoms that are subject to the forces from the
density-guided simulation and contribute to the simulated density.


density-guided-simulation-similarity-measure

(inner-product) Similarity measure between the density that is
calculated from the atom positions and the reference density.

inner-product

Takes the sum of the product of reference density and simulated density
voxel values.


relative-entropy

Uses the negative relative entropy (or Kullback-Leibler divergence)
between reference density and simulated density as similarity measure.
Negative density values are ignored.


cross-correlation

Uses the Pearson correlation coefficient between reference density and
simulated density as similarity measure.


density-guided-simulation-atom-spreading-weight

(unity) Determines the multiplication factor for the Gaussian kernel
when spreading atoms on the grid.

unity

Every atom in the density fitting group is assigned the same unit
factor.


mass

Atoms contribute to the simulated density proportional to their mass.


charge

Atoms contribute to the simulated density proportional to their charge.


density-guided-simulation-force-constant

(1e+09) \[kJ mol^{-1}\] The scaling factor for density-guided
simulation forces. May also be negative.


density-guided-simulation-gaussian-transform-spreading-width

(0.2) \[nm\] The Gaussian RMS width for the spread kernel for the
simulated density.


density-guided-simulation-gaussian-transform-spreading-range-in-multiples-of-width

\(4\) The range after which the gaussian is cut off in multiples of the
Gaussian RMS width described above.


density-guided-simulation-reference-density-filename

(reference.mrc) Reference density file name using an absolute path or a
path relative to the to the folder from which `gmx mdrun` is called.


density-guided-simulation-nst

\(1\) Interval in steps at which the density fitting forces are
evaluated and applied. The forces are scaled by this number when applied
(See the [reference manual]() for details).


density-guided-simulation-normalize-densities

(true) Normalize the sum of density voxel values to one for the
reference density as well as the simulated density.


density-guided-simulation-adaptive-force-scaling

(false) Adapt the force constant to ensure a steady increase in
similarity between simulated and reference density.

> Do not use adaptive force scaling.

true

Use adaptive force scaling.


density-guided-simulation-adaptive-force-scaling-time-constant

\(4\) \[ps\] Couple force constant to increase in similarity with
reference density with this time constant. Larger times result in looser
coupling.


density-guided-simulation-shift-vector

(0,0,0) \[nm\] Add this vector to all atoms in the
density-guided-simulation-group before calculating forces and energies
for density-guided simulations. Affects only the density-guided
simulation forces and energies. Corresponds to a shift of the input
density in the opposite direction by (-1) \*
density-guided-simulation-shift-vector.


density-guided-simulation-transformation-matrix

(1,0,0,0,1,0,0,0,1) Multiply all atoms with this matrix in the
density-guided-simulation-group before calculating forces and energies
for density-guided simulations. Affects only the density-guided
simulation forces and energies. Corresponds to a transformation of the
input density by the inverse of this matrix. The matrix is given in
row-major order. This option allows, e.g., rotation of the
density-guided atom group around the z-axis by $\theta$ degrees by
using the following input:
$(\cos \theta , -\sin \theta , 0 , \sin \theta , \cos \theta , 0 , 0 , 0 , 1)$
.


### QM/MM simulations with CP2K Interface

These options enable and control the calculation and application of
additional QM/MM forces that are computed by the CP2K package if it is
linked into GROMACS. For further details about
QM/MM interface implementation follow `qmmm`.

qmmm-cp2k-active

(false) Activate QM/MM simulations. Requires CP2K to be linked with
GROMACS


qmmm-cp2k-qmgroup

(System) Index group with atoms that are treated with QM.


qmmm-cp2k-qmmethod

(PBE) Method used to describe the QM part of the system.

PBE

DFT using PBE functional and DZVP-MOLOPT basis set.


PBE-D3

DFT using PBE functional, DZVP-MOLOPT basis set and D3 dispersion
correction.


BLYP

DFT using BLYP functional and DZVP-MOLOPT basis set.


BLYP-D3

DFT using BLYP functional, DZVP-MOLOPT basis set and D3 dispersion
correction.


B3LYP

DFT using B3LYP hybrid functional and 6-31G\* basis set.


B3LYP-D3

DFT using B3LYP hybrid functional, 6-31G\* basis set and D3 dispersion
correction.


PBE0

DFT using PBE0 hybrid functional and 6-31G\* basis set.


PBE0-D3

DFT using PBE0 hybrid functional, 6-31G\* basis set and D3 dispersion
correction.


CAM-B3LYP

DFT using CAM-B3LYP range-separated functional and 6-31G\* basis set.


CAM-B3LYP-D3

DFT using CAM-B3LYP range-separated functional, 6-31G\* basis set and D3
dispersion correction.


WB97X

DFT using wB97X range-separated functional and 6-31G\* basis set.


WB97X-D3

DFT using wB97X-D3 range-separated functional with built-in D3
dispersion correction and 6-31G\* basis set.


INPUT

Provide an external input file for CP2K when running `gmx grompp` with
the `-qmi` command-line option. External input files are subject to the
limitations that are described in `qmmm`.


qmmm-cp2k-qmcharge

0)  Total charge of the QM part.


qmmm-cp2k-qmmultiplicity

1)  Multiplicity or spin-state of QM part. Default value 1 means singlet
    state.


qmmm-cp2k-qmfilenames

() Names of the CP2K files that will be generated during the simulation.
When using the default, empty, value the name of the simulation input
file will be used with an additional `_cp2k` suffix.


### Collective variables (Colvars) module

These options enable and control the features provided by the collective
variables ([Colvars](https://colvars.github.io/)) module, a software
library for enhanced sampling methods in molecular simulations. The
Colvars module is described in ref. `195 <refFiorin13>` as well as other
references that are reported in the log file when the corresponding
features are used. For further details about Colvars interface
implementation follow `colvars`.

colvars-active

(false) Activate Colvars computation in the current run. Requires that
the Colvars library was compiled with GROMACS,
which is the default in a typical installation.


colvars-configfile

Name of the Colvars configuration file, using options specific to
Colvars that are documented at:
<https://colvars.github.io/gromacs-2026/colvars-refman-gromacs.html>.
The file name can be either an absolute path, or a path relative to the
working directory when `gmx grompp` is called.


colvars-seed

(-1) \[integer\] Seed used to initialize the random generator associated
with certain stochastic methods implemented within Colvars. The default
value of -1 generates a random seed.


The current implementation of the
Colvars-GROMACS interface gathers the relevant
atomic coordinates on one MPI rank, where all collective variables and
their forces are computed. Take this fact into account when choosing how
many atoms to include in selections.

### NNP/MM simulations with neural network potentials in the NNPot interface

These options enable and control the calculation and application of
additional forces derived from neural network potentials via the nnpot
interface, if GROMACS is built with LibTorch
support. For further details about NNP/MM interface implementation
follow `nnpot`.

nnpot-active

(false) Activate NNP/MM simulations via the NNPot interface. Requires
GROMACS to be built with LibTorch support.


nnpot-modelfile

(model.pt) Path to a TorchScript-compiled model, either absolute or
relative to the simulation directory.


nnpot-input-group

(System) Index group defining the input atoms for the NNP subsystem.
Defaults to `System`, which performs a pure NNP simulation.


nnpot-embedding

(mechanical) Type of embedding scheme to use for NNP/MM simulations.

mechanical

Mechanical embedding scheme, NNP-MM interactions are treated
classically.


electrostatic-model

Electrostatic embedding scheme, NNP-MM interactions are computed by the
NNP model. In this case, the NNP model is expected to return the total
energy, as well as forces on NNP and MM atoms.


nnpot-model-input\[1-9\]

Names of input fields to be filled in by
GROMACS at each step and passed to the NNP
model. Up to 9 inputs can be specified:

atom-positions

Positions of the NNP atoms specified by `nnpot-input-group`.


atom-numbers

Atomic numbers of the NNP atoms specified by `nnpot-input-group`.


atom-pairs

Pairs of NNP atoms specified by `nnpot-input-group`, filter by the
`pair-cutoff`, as a vector of shape (N_pairs, 2).


pair-shifts

Periodic box shift vectors for the NNP atom pairs, as a vector of shape
(N_pairs, 3).


atom-positions-mm

Positions of the MM atoms surrounding the NNP region, needed for
`electrostatic-model` embedding.


atom-charges-mm

Charges of the MM atoms surrounding the NNP region, needed for
`electrostatic-model` embedding.


nnp-charge

Total charge of the NNP atoms.


box

Simulation box vectors.


pbc

Periodic boundary conditions.


pair-cutoff

(0.0) \[nm\] Cutoff distance for pairs of NNP atoms. Positive cutoff
value must be specified when requesting `atom-pairs` input.


nnpot-link-type

\(H\) Type of link atoms to be used in the NNP/MM simulation, specified
by element symbol. Defaults to hydrogen (H).


nnpot-link-distance

(0.1) \[nm\] Distance between link atom and the bonded MM atom.


### Fast Multipole Method (FMM) Interface

These options enable an FMM backend and control its settings. For
further details about the FMM interface implementation follow `fmm`.

fmm-backend

(inactive) Activate an FMM library for electrostatics.

inactive

No FMM library is used.


exafmm

Use ExaFMM as FMM library.


fmsolvr

Use FMSolvr as FMM library.


#### Options controlling ExaFmm

fmm-exafmm-direct-provider

(GROMACS) Selects direct interaction provider.

GROMACS

Use GROMACS to compute direct interactions.


FMM

Use ExaFMM to compute direct interactions.


fmm-exafmm-direct-range

\(2\) Specifies whether ExaFMM includes only immediate or also
second-nearest neighboring cells in direct interactions. Accepts `1` or
`2`. Must be `2` when using GROMACS as direct
provider.


fmm-exafmm-order

6)  Multipole expansion order. Must be greater than `0`.


fmm-exafmm-tree-type

(uniform) Tree structure used by ExaFMM.

uniform

Use uniform tree structure.


adaptive

Use adaptive tree structure.


fmm-exafmm-tree-depth

\(0\) Tree depth for uniform tree. Required when using `FMM` as direct
provider. Must not be set for adaptive tree.


fmm-exafmm-max-particles-per-cell

\(0\) Maximum particles per cell for adaptive tree. Required and
positive when `fmm-exafmm-tree-type = adaptive`. Must not be set for
uniform tree.


#### Options controlling FMSolvr

fmm-fmsolvr-order

8)  Multipole expansion order. Must be greater than `0`.


fmm-fmsolvr-tree-depth

3)  Tree depth controlling spatial subdivision.


fmm-fmsolvr-direct-provider

(FMM) Selects direct interaction provider.

GROMACS

Use GROMACS to compute direct interactions.


FMM

Use FMSolvr to compute direct interactions.


fmm-fmsolvr-direct-range

1)  Direct interaction settings. Must be `1`, as FMSolvr currently
    supports direct interactions only with immediate neighboring cells.


fmm-fmsolvr-dipole-compensation

(yes) Enables dipole compensation.


fmm-fmsolvr-sparse

(no) Enable performance optimizations for sparse simulation systems with
a lot of vacuum.


### User defined thingies

user1-grps


user2-grps


userint1 (0)


userint2 (0)


userint3 (0)


userint4 (0)


userreal1 (0)


userreal2 (0)


userreal3 (0)


userreal4 (0)

These you can use if you modify code. You can pass integers and reals
and groups to your subroutine. Check the inputrec definition in
`src/gromacs/mdtypes/inputrec.h`


### Removed features

These features have been removed from GROMACS,
but so that old `mdp` and `tpr` files cannot be mistakenly misused, we
still parse these options. `gmx grompp` and `gmx mdrun` will issue a
fatal error if they are set.

adress

(no)


implicit-solvent

(no)


