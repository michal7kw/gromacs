# Replica exchange

Replica exchange molecular dynamics (REMD) is a method that can be used
to speed up the sampling of any type of simulation, especially if
conformations are separated by relatively high energy barriers. It
involves simulating multiple replicas of the same system at different
temperatures and randomly exchanging the complete state of two replicas
at regular intervals with the probability:

$$
P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
\left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2)
\right] \right)
$$

where $T_1$ and $T_2$ are the reference temperatures and $U_1$ and
$U_2$ are the instantaneous potential energies of replicas 1 and 2
respectively. After exchange the velocities are scaled by
$(T_1/T_2)^{\pm0.5}$ and a neighbor search is performed the next step.
This combines the fast sampling and frequent barrier-crossing of the
highest temperature with correct Boltzmann sampling at all the different
temperatures `60 <refHukushima96a>`, `61 <refSugita99>`. We only attempt
exchanges for neighboring temperatures as the probability decreases very
rapidly with the temperature difference. One should not attempt
exchanges for all possible pairs in one step. If, for instance, replicas
1 and 2 would exchange, the chance of exchange for replicas 2 and 3 not
only depends on the energies of replicas 2 and 3, but also on the energy
of replica 1. In GROMACS this is solved by
attempting exchange for all *odd* pairs on *odd* attempts and for all
*even* pairs on *even* attempts. If we have four replicas: 0, 1, 2 and
3, ordered in temperature and we attempt exchange every 1000 steps,
pairs 0-1 and 2-3 will be tried at steps 1000, 3000 etc. and pair 1-2 at
steps 2000, 4000 etc.

How should one choose the temperatures? The energy difference can be
written as:

$$
U_1 - U_2 =  N_{df} \frac{c}{2} k_B (T_1 - T_2)
$$

where $N_{df}$ is the total number of degrees of freedom of one
replica and $c$ is 1 for harmonic potentials and around 2 for
protein/water systems. If $T_2 = (1+\epsilon) T_1$ the probability
becomes:

$$
P(1 \leftrightarrow 2)
= \exp\left( -\frac{\epsilon^2 c\,N_{df}}{2 (1+\epsilon)} \right)
\approx \exp\left(-\epsilon^2 \frac{c}{2} N_{df} \right)
$$

Thus for a probability of $e^{-2}\approx 0.135$ one obtains
$\epsilon \approx 2/\sqrt{c\,N_{df}}$. With all bonds constrained one
has $N_{df} \approx 2\, N_{atoms}$ and thus for $c$ = 2 one should
choose $\epsilon$ as $1/\sqrt{N_{atoms}}$. However there is one
problem when using pressure coupling. The density at higher temperatures
will decrease, leading to higher energy `62 <refSeibert2005a>`, which
should be taken into account. Using a so-called [REMD
calculator](https://virtualchemistry.org/remd-temperature-generator/),
you can type in the temperature range and the number of atoms. The tool
then proposes a set of temperatures.

An extension to the REMD for the isobaric-isothermal ensemble was
proposed by Okabe et al. `63 <refOkabe2001a>`. In this work the exchange
probability is modified to:

$$
P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
\left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2) +
\left(\frac{P_1}{k_B T_1} - \frac{P_2}{k_B T_2}\right)\left(V_1-V_2\right)
\right] \right)
$$

where $P_1$ and $P_2$ are the respective reference pressures and
$V_1$ and $V_2$ are the respective instantaneous volumes in the
simulations. In most cases the differences in volume are so small that
the second term is negligible. It only plays a role when the difference
between $P_1$ and $P_2$ is large or in phase transitions.

Hamiltonian replica exchange is also supported in
GROMACS. In Hamiltonian replica exchange, each
replica has a different Hamiltonian, defined by the free energy pathway
specified for the simulation. The exchange probability to maintain the
correct ensemble probabilities is:

$$
P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
\frac{1}{k_B T} (U_1(x_1) - U_1(x_2) + U_2(x_2) - U_2(x_1))
\right]\right)
$$

The separate Hamiltonians are defined by the free energy functionality
of GROMACS, with swaps made between the
different values of $\lambda$ defined in the mdp file.

Hamiltonian and temperature replica exchange can also be performed
simultaneously `64 <refChodera2011>`, using the acceptance criteria:

$$
P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
\frac{U_1(x_1) - U_1(x_2)}{k_B T_1} + \frac{U_2(x_2) - U_2(x_1)}{k_B T_2}
\right] \right)
$$

Gibbs sampling replica exchange has also been implemented in
GROMACS `64 <refChodera2011>`. In Gibbs
sampling replica exchange, all possible pairs are tested for exchange,
allowing swaps between replicas that are not neighbors.

Gibbs sampling replica exchange requires no additional potential energy
calculations. However there is an additional communication cost in Gibbs
sampling replica exchange, as for some permutations, more than one round
of swaps must take place. In some cases, this extra communication cost
might affect the efficiency.

All replica exchange variants are set using `mdp` options and performed
using the `mdrun ` program. It will only work when MPI is
installed, due to the inherent parallelism in the algorithm. For
efficiency each replica can run on a separate rank. See the manual page
of `mdrun ` on how to use these multinode features.
