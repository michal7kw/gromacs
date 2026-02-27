# Free energy interactions

This section describes the $\lambda$-dependence of the potentials used
for free energy calculations (see sec. `fecalc`). All common types of
potentials and constraints can be interpolated smoothly from state A
($\lambda=0$) to state B ($\lambda=1$) and vice versa. All bonded
interactions are interpolated by linear interpolation of the interaction
parameters. Non-bonded interactions can be interpolated linearly or via
soft-core interactions.

Starting in GROMACS 4.6, $\lambda$ is a
vector, allowing different components of the free energy transformation
to be carried out at different rates. Coulomb, Lennard-Jones, bonded,
and restraint terms can all be controlled independently, as described in
the `mdp` options.

## Harmonic potentials

The example given here is for the bond potential, which is harmonic in
GROMACS. However, these equations apply to the
angle potential and the improper dihedral potential as well.

$$
\begin{aligned}
\begin{aligned}
V_b     &=&{\frac{1}{2}}\left[{(1-{\lambda})}k_b^A +
{\lambda}k_b^B\right] \left[b - {(1-{\lambda})}b_0^A - {\lambda}b_0^B\right]^2  \\
{\frac{\partial V_b}{\partial {\lambda}}}&=&{\frac{1}{2}}(k_b^B-k_b^A)
\left[b - {(1-{\lambda})}b_0^A + {\lambda}b_0^B\right]^2 +
\nonumber\\
& & \phantom{{\frac{1}{2}}}(b_0^A-b_0^B) \left[b - {(1-{\lambda})}b_0^A -{\lambda}b_0^B\right]
\left[{(1-{\lambda})}k_b^A + {\lambda}k_b^B \right]\end{aligned}
\end{aligned}
$$

## GROMOS-96 bonds and angles

Fourth-power bond stretching and cosine-based angle potentials are
interpolated by linear interpolation of the force constant and the
equilibrium position. Formulas are not given here.

## Proper dihedrals

For the proper dihedrals, the equations are somewhat more complicated:

$$
\begin{aligned}
\begin{aligned}
V_d     &=&\left[{(1-{\lambda})}k_d^A + {\lambda}k_d^B \right]
\left( 1+ \cos\left[n_{\phi} \phi -
{(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
\right]\right)\\
{\frac{\partial V_d}{\partial {\lambda}}}&=&(k_d^B-k_d^A)
\left( 1+ \cos
\left[
n_{\phi} \phi- {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
\right]
\right) +
\nonumber\\
&&(\phi_s^B - \phi_s^A) \left[{(1-{\lambda})}k_d^A - {\lambda}k_d^B\right]
\sin\left[  n_{\phi}\phi - {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B \right]\end{aligned}
\end{aligned}
$$

**Note:** that the multiplicity $n_{\phi}$ can not be parameterized
because the function should remain periodic on the interval
$[0,2\pi]$.

## Tabulated bonded interactions

For tabulated bonded interactions only the force constant can
interpolated:

$$
\begin{aligned}
\begin{aligned}
V  &=& ({(1-{\lambda})}k^A + {\lambda}k^B) \, f \\
{\frac{\partial V}{\partial {\lambda}}} &=& (k^B - k^A) \, f\end{aligned}
\end{aligned}
$$

## Coulomb interaction

The Coulomb interaction between two particles of which the charge varies
with ${\lambda}$ is:

$$
\begin{aligned}
\begin{aligned}
V_c &=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
{\frac{\partial V_c}{\partial {\lambda}}}&=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[- q_i^A q_j^A + q_i^B q_j^B\right]\end{aligned}
\end{aligned}
$$

where $f = \frac{1}{4\pi \varepsilon_0} = {138.935\,458}$ (see
chapter `defunits`).

## Coulomb interaction with reaction field

The Coulomb interaction including a reaction field, between two
particles of which the charge varies with ${\lambda}$ is:

$$
\begin{aligned}
\begin{aligned}
V_c     &=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
\left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
{\frac{\partial V_c}{\partial {\lambda}}}&=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
\left[- q_i^A q_j^A + q_i^B q_j^B\right]
\end{aligned}
\end{aligned}
$$

**Note** that the constants $k_{rf}$ and $c_{rf}$ are defined using
the dielectric constant ${\varepsilon_{rf}}$ of the medium (see
sec. `coulrf`).

## Lennard-Jones interaction

For the Lennard-Jones interaction between two particles of which the
*atom type* varies with ${\lambda}$ we can write:

$$
\begin{aligned}
\begin{aligned}
V_{LJ}  &=&     \frac{{(1-{\lambda})}C_{12}^A + {\lambda}\, C_{12}^B}{{r_{ij}}^{12}} -
\frac{{(1-{\lambda})}C_6^A + {\lambda}\, C_6^B}{{r_{ij}}^6}   \\
{\frac{\partial V_{LJ}}{\partial {\lambda}}}&=&\frac{C_{12}^B - C_{12}^A}{{r_{ij}}^{12}} -
\frac{C_6^B - C_6^A}{{r_{ij}}^6}
\end{aligned}
\end{aligned}
$$

It should be noted that it is also possible to express a pathway from
state A to state B using $\sigma$ and $\epsilon$ (see
`eqn. %s <eqnsigeps>`). It may seem to make sense physically to vary the
force field parameters $\sigma$ and $\epsilon$ rather than the
derived parameters $C_{12}$ and $C_{6}$. However, the difference
between the pathways in parameter space is not large, and the free
energy itself does not depend on the pathway, so we use the simple
formulation presented above.

## Kinetic Energy

When the mass of a particle changes, there is also a contribution of the
kinetic energy to the free energy (note that we can not write the
momentum $\mathbf{p}$ as $m \mathbf{v}$, since that would result in
the sign of ${\frac{\partial E_k}{\partial {\lambda}}}$ being
incorrect `99 <refGunsteren98a>`):

$$
\begin{aligned}
\begin{aligned}
E_k      &=&     {\frac{1}{2}}\frac{\mathbf{p}^2}{{(1-{\lambda})}m^A + {\lambda}m^B}        \\
{\frac{\partial E_k}{\partial {\lambda}}}&=&    -{\frac{1}{2}}\frac{\mathbf{p}^2(m^B-m^A)}{({(1-{\lambda})}m^A + {\lambda}m^B)^2}\end{aligned}
\end{aligned}
$$

after taking the derivative, we *can* insert $\mathbf{p}$ =
$m \mathbf{v}$, such that:

$$
{\frac{\partial E_k}{\partial {\lambda}}}~=~    -{\frac{1}{2}}\mathbf{v}^2(m^B-m^A)
$$

## Constraints

The constraints are formally part of the Hamiltonian, and therefore they
give a contribution to the free energy. In
GROMACS this can be calculated using the LINCS
or the SHAKE algorithm. If we have $k = 1 \ldots K$ constraint
equations $g_k$ for LINCS, then

$$
g_k     =       | \mathbf{r}_{k} | - d_{k}
$$

where $\mathbf{r}_k$ is the displacement vector between two particles
and $d_k$ is the constraint distance between the two particles. We can
express the fact that the constraint distance has a ${\lambda}$
dependency by

$$
d_k     =       {(1-{\lambda})}d_{k}^A + {\lambda}d_k^B
$$

Thus the ${\lambda}$-dependent constraint equation is

$$
g_k     =       | \mathbf{r}_{k} | - \left({(1-{\lambda})}d_{k}^A + {\lambda}d_k^B\right).
$$

The (zero) contribution $G$ to the Hamiltonian from the constraints
(using Lagrange multipliers $\lambda_k$, which are logically distinct
from the free-energy ${\lambda}$) is

$$
\begin{aligned}
\begin{aligned}
G           &=&     \sum^K_k \lambda_k g_k    \\
{\frac{\partial G}{\partial {\lambda}}}    &=&     \frac{\partial G}{\partial d_k} {\frac{\partial d_k}{\partial {\lambda}}} \\
&=&     - \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}
\end{aligned}
$$

For SHAKE, the constraint equations are

$$
g_k     =       \mathbf{r}_{k}^2 - d_{k}^2
$$

with $d_k$ as before, so

$$
\begin{aligned}
{\frac{\partial G}{\partial {\lambda}}}    &=&     -2 \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}
$$

## Soft-core interactions: Beutler *et al.*

![plots/softcore.*](plots/softcore.*)

*Soft-core interactions at $λ = 0.5$, with $p = 2$ and $C6A = C12A = C6B = C12B = 1$.*

In a free-energy calculation where particles grow out of nothing, or
particles disappear, using the simple linear interpolation of the
Lennard-Jones and Coulomb potentials as described in
`Equations %s <eqdVljdlambda>` and `%s <eqdVcoulombdlambda>` may lead to
poor convergence. When the particles have nearly disappeared, or are
close to appearing (at ${\lambda}$ close to 0 or 1), the interaction
energy will be weak enough for particles to get very close to each
other, leading to large fluctuations in the measured values of
$\partial V/\partial {\lambda}$ (which, because of the simple linear
interpolation, depends on the potentials at both the endpoints of
${\lambda}$).

To circumvent these problems, the singularities in the potentials need
to be removed. This can be done by modifying the regular Lennard-Jones
and Coulomb potentials with “soft-core” potentials that limit the
energies and forces involved at ${\lambda}$ values between 0 and 1,
but not *at* ${\lambda}=0$ or 1.

In GROMACS the soft-core potentials $V_{sc}$
are shifted versions of the regular potentials, so that the singularity
in the potential and its derivatives at $r=0$ is never reached. This
formulation was introduced by Beutler *et al.*`100 <refBeutler94>`:

$$
\begin{aligned}
\begin{aligned}
V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
\\
r_A &=& \left(\alpha \sigma_A^6 {\lambda}^p + r^6 \right)^\frac{1}{6}
\\
r_B &=& \left(\alpha \sigma_B^6 {(1-{\lambda})}^p + r^6 \right)^\frac{1}{6}\end{aligned}
\end{aligned}
$$

where $V^A$ and $V^B$ are the normal “hard core” Van der Waals or
electrostatic potentials in state A (${\lambda}=0$) and state B
(${\lambda}=1$) respectively, $\alpha$ is the soft-core parameter
(set with `sc_alpha` in the `mdp` file), $p$ is the soft-core
${\lambda}$ power (set with `sc_power`), $\sigma$ is the radius of
the interaction, which is $(C_{12}/C_6)^{1/6}$ or an input parameter
(`sc_sigma`) when $C_6$ or $C_{12}$ is zero. Beutler *et
al.*`100 <refBeutler94>` probed various combinations of the $r$ power
values for the Lennard-Jones and Coulombic interactions.
GROMACS uses $r^6$ for both, van der Waals
and electrostatic interactions.

For intermediate ${\lambda}$, $r_A$ and $r_B$ alter the
interactions very little for $r > \alpha^{1/6} \sigma$ and quickly
switch the soft-core interaction to an almost constant value for smaller
$r$ (`Fig. %s <fig-softcore>`). The force is:

$$
F_{sc}(r) = -\frac{\partial V_{sc}(r)}{\partial r} =
{(1-{\lambda})}F^A(r_A) \left(\frac{r}{r_A}\right)^5 +
{\lambda}F^B(r_B) \left(\frac{r}{r_B}\right)^5
$$

where $F^A$ and $F^B$ are the “hard core” forces. The contribution
to the derivative of the free energy is:

$$
\begin{aligned}
\begin{aligned}
{\frac{\partial V_{sc}(r)}{\partial {\lambda}}} & = &
V^B(r_B) -V^A(r_A)  +
{(1-{\lambda})}\frac{\partial V^A(r_A)}{\partial r_A}
\frac{\partial r_A}{\partial {\lambda}} +
{\lambda}\frac{\partial V^B(r_B)}{\partial r_B}
\frac{\partial r_B}{\partial {\lambda}}
\nonumber\\
&=&
V^B(r_B) -V^A(r_A)  + \nonumber \\
& &
\frac{p \alpha}{6}
\left[ {\lambda}F^B(r_B) r^{-5}_B \sigma_B^6 {(1-{\lambda})}^{p-1} -
{(1-{\lambda})}F^A(r_A) r^{-5}_A \sigma_A^6 {\lambda}^{p-1} \right]\end{aligned}
\end{aligned}
$$

The original GROMOS Lennard-Jones soft-core function
`100 <refBeutler94>` uses $p=2$, but $p=1$ gives a smoother
$\partial H/\partial{\lambda}$ curve. Another issue that should be
considered is the soft-core effect of hydrogens without Lennard-Jones
interaction. Their soft-core $\sigma$ is set with `sc_sigma`. These
hydrogens produce peaks in $\partial H/\partial{\lambda}$ at
${\lambda}$ is 0 and/or 1 for $p=1$ and close to 0 and/or 1 with
$p=2$. Lowering `sc_sigma` will decrease this effect, but it will also
increase the interactions with hydrogens relative to the other
interactions in the soft-core state.

When soft-core potentials are selected (by setting `sc_alpha >0`), and
the Coulomb and Lennard-Jones potentials are turned on or off
sequentially, then the Coulombic interaction is turned off linearly,
rather than using soft-core interactions, which should be less
statistically noisy in most cases. This behavior can be overwritten by
setting `sc-coul=yes`. Note that `sc-coul` is only taken into account
when lambda states are used, and you can still turn off soft-core
interactions by setting `sc-alpha=0`. Additionally, the soft-core
interaction potential is only applied when either the A or B state has
zero interaction potential. If both A and B states have nonzero
interaction potential, default linear scaling described above is used.
When both Coulombic and Lennard-Jones interactions are turned off
simultaneously, a soft-core potential is used, and a hydrogen is being
introduced or deleted, the sigma is set to `sc-sigma-min`, which itself
defaults to `sc-sigma-default`.

## Soft-core interactions: Gapsys *et al.*

In this section we describe the functional form and parameters for the
soft-cored non-bonded interactions using the formalism by Gapsys *et
al.* `183 <refGapsys2012>`.

The Gapsys *et al.* soft-core is formulated to act on the level of van
der Waals and electrostatic forces: the non-bonded interactions are
linearized at a point defined as, $r_{scLJ}$ or $r_{scQ}$,
respectively. The linearization point depends on the state of the system
as controlled by the $\lambda$ parameter and two parameters
$\alpha_Q$ (set with `sc-gapsys-scale-linpoint-q`) and $\alpha_{LJ}$
(set with `sc-gapsys-scale-linpoint-lj`). The dependence on $\lambda$
guarantees that the end-states are properly represented by their
hard-core potentials. `Fig. %s <fig-gapsyssc>` illustrates the behaviour
of the linearization point, forces and integrated potential energies
with respect to the parameters $\alpha_Q$ and $\alpha_{LJ}$. The
optimal choices of the parameter values have been systematically
explored in `183 <refGapsys2012>`. These recommended values are set by
default when `sc-function=gapsys` is selected:
`sc-gapsys-scale-linpoint-q=0.3` and `sc-gapsys-scale-linpoint-lj=0.85`.

![plots/gapsys-sc.*](plots/gapsys-sc.*)

*Illustration of the soft-core parameter influence on the linearization point (top row), forces (middle row) and energies (bottom row) for van der Waals (left column) and electrostatic interactions (right column). The case of two interacting atoms is considered. In state A both atoms have charges of 0.5 and $σ = 0.3$ nm, $ϵ = 0.5$ kJ/mol. In state B all the non-bonded interactions are set to zero. The parameter $λ$ is set to 0.5 and electrostatic interaction cutoff is 1 nm.*

The parameter $\alpha_{LJ}$ is a unitless scaling factor in the range
$[0,1)$. It scales the position of the point from which the van der
Waals force will be linearized. The linearization of the force is
allowed in the range $[0,F_{min}^{LJ})$, where setting
$\alpha_{LJ}=0$ results in a standard hard-core van der Waals
interaction. Setting $\alpha_{LJ}$ closer to 1 brings the force
linearization point towards the minimum in the Lennard-Jones force curve
($F_{min}^{LJ}$). This construct allows retaining the repulsion
between two particles with non-zero C12 parameter at any $\lambda$
value.

The parameter $\alpha_{Q}$ has a unit of $\frac{nm}{e^2}$ and is
defined in the range $[0,\inf)$. It scales the position of the point
from which the Coulombic force will be linearized. Even though in theory
$\alpha_{Q}$ can be set to an arbitrarily large value, algorithmically
the linearization point for the force is bound in the range
$[0,F_{rcoul}^{Q})$, where setting $\alpha_{Q}=0$ results in a
standard hard-core Coulombic interaction. Setting $\alpha_{Q}$ to a
larger value softens the Coulombic force.

In all the notations below, for simplicity, the distance between two
atoms $i$ and $j$ is noted as $r$, i.e. $r=r_{ij}$.

### Forces: van der Waals interactions

$$
\begin{aligned}
\begin{aligned}
\mathbf{F}_{ij}^{LJ}(\mathbf{r})=\begin{cases}
(\frac{12C_{ij}^{(12)}}{r^{13}} - \frac{6C_{ij}^{(6)}}{r^7})\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$}
\\
\frac{d\mathbf{F}_{ij}^{LJ}}{dr}_{r=r_{scLJ}}r + \mathbf{F}_{ij}^{LJ}(r_{scLJ}), & \mbox{if } \mbox{ $r<r_{scLJ}$}
\end{cases}\end{aligned}
\end{aligned}
$$

where the switching point between the soft and hard-core Lennard-Jones
forces
$r_{scLJ} = \alpha_{LJ}(\frac{26}{7}\sigma^6\lambda)^{\frac{1}{6}}$
for state A, and
$r_{scLJ} = \alpha_{LJ}(\frac{26}{7}\sigma^6(1-\lambda))^{\frac{1}{6}}$
for state B. In analogy to the Beutler *et al.* soft core version,
$\sigma$ is the radius of the interaction, which is
$(C_{12}/C_6)^{1/6}$ or an input parameter (set with
`sc-sigma-LJ-gapsys`) when C6 or C12 is zero. The default value for this
parameter is `sc-sigma-LJ-gapsys=0.3`.

Explicit expression:

$$
\begin{aligned}
\begin{aligned}
\mathbf{F}_{LJ}(\mathbf{r})=\begin{cases}
\left(\frac{12C^{(12)}}{r^{13}} - \frac{6C^{(6)}}{r^7}\right)\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$}
\\
\left(-\frac{156C^{(12)}}{r_{scLJ}^{14}} + \frac{42C^{(6)}}{r_{scLJ}^{8}}\right)\mathbf{r} + \frac{168C^{(12)}}{r_{scLJ}^{13}} - \frac{48C^{(6)}}{r_{scLJ}^{7}}, & \mbox{if } \mbox{ $r<r_{scLJ}$}
\end{cases}\end{aligned}
\end{aligned}
$$

### Forces: Coulomb interactions

$$
\begin{aligned}
\begin{aligned}
\mathbf{F}_{ij}^{Q}(\mathbf{r})=\begin{cases}
\frac{q_{i}q_{j}}{4{\pi}{\varepsilon_0}{\varepsilon_r}r^{2}}\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$}
\\
\frac{d\mathbf{F}_{ij}^{Q}}{dr}_{r=r_{scQ}}r + \mathbf{F}_{ij}^{Q}(r_{scQ}), & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
\\
\frac{d\mathbf{F}_{ij}^{Q}}{dr}_{r=r_{cutoffQ}}r + \mathbf{F}_{ij}^{Q}(r_{cutoffQ}), & \mbox{if } \mbox{ $r < r^{scQ} \geq r_{cutoffQ}$}
\end{cases}\end{aligned}
\end{aligned}
$$

where the switching point $r^{sc}$ between the soft and hard-core
electrostatic forces is
$r_{scQ} = \alpha_Q(1+|q_iq_j|)\lambda^{\frac{1}{6}}$ for state A, and
$r_{scQ} = \alpha_Q(1+|q_iq_j|)(1-\lambda)^{\frac{1}{6}}$ for state B.
The $\lambda$ dependence of the linearization point for both van der
Waals and Coulombic interactions is of the same power $1/6$.

Explicit expression:

$$
\begin{aligned}
\begin{aligned}
\mathbf{F}_{Q}(\mathbf{r})=\begin{cases}
\frac{q_iq_j}{4{\pi}{\varepsilon_0}{\varepsilon_r}r^{2}}\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$}
\\
\frac{1}{4{\pi}{\varepsilon_0}{\varepsilon_r}}\big( -\frac{2q_{i}q_{j}}{r_{sc}^3}\mathbf{r} + \frac{3q_iq_j}{r_{sc}^2} \big), & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
\\
\frac{1}{4{\pi}{\varepsilon_0}{\varepsilon_r}}\big( -\frac{2q_{i}q_{j}}{r_{cutoffQ}^3}\mathbf{r} + \frac{3q_iq_j}{r_{cutoffQ}^2} \big), & \mbox{if } \mbox{ $r < r_{scQ} \geq r_{cutoffQ}$}                        \end{cases}\end{aligned}
\end{aligned}
$$

### Energies: van der Waals interactions

Explicit definition of energies:

$$
\begin{aligned}
\begin{aligned}
V_{LJ}(r)=\begin{cases}
\frac{C^{(12)}}{r^{12}} - \frac{C^{(6)}}{r^6}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$}
\\
\left(\frac{78C^{(12)}}{r_{scLJ}^{14}} - \frac{21C^{(6)}}{r_{scLJ}^{8}}\right)r^2 - \left(\frac{168C^{(12)}}{r_{scLJ}^{13}} - \frac{48C^{(6)}}{r_{scLJ}^{7}}\right)r
+ \frac{91C^{(12)}}{r_{scLJ}^{12}} - \frac{28C^{(6)}}{r_{scLJ}^{6}}, & \mbox{if } \mbox{ $r<r_{scLJ}$}
\end{cases}\end{aligned}
\end{aligned}
$$

### Energies: Coulomb interactions

$$
\begin{aligned}
\begin{aligned}
V_{Q}(r)=\begin{cases}
\frac{q_{i}q_{j}}{4{\pi}{\varepsilon_0}{\varepsilon_r}r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$}
\\
\frac{q_{i}q_{j}}{r_{scQ}^3}r^2 - \frac{3q_iq_j}{r_{scQ}^2}r + \frac{3q_iq_j}{r_{scQ}}, & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
\\
\frac{q_{i}q_{j}}{r_{cutoffQ}^3}r^2 - \frac{3q_iq_j}{r_{cutoffQ}^2}r + \frac{3q_iq_j}{r_{cutoffQ}}, & \mbox{if } \mbox{ $r < r_{scQ} \geq r_{cutoffQ}$}
\end{cases}\end{aligned}
\end{aligned}
$$

### $\partial H / \partial \lambda$: van der Waals interactions

Here we provide the explicit expressions of
$\partial H/ \partial \lambda$ for Lennard-Jones potential, when
$r<r_{scLJ}$. For simplicity, in the expression below we use the
notation $r_{scLJ_A}=r_{scA}$ and $r_{scLJ_B}=r_{scB}$.

$$
\begin{aligned}
\begin{aligned}
\frac{\partial{H}}{\partial{\lambda}} &= V_{LJ}^B(r) - V_{LJ}^A(r) + (1-\lambda)\frac{\partial{V_{LJ}^A(r)}}{\partial{\lambda}} + \lambda\frac{\partial{V_{LJ}^B(r)}}{\partial{\lambda}} \\
& =  \left(\frac{78C^{(12)}_B}{r_{scB}^{14}} - \frac{21C^{(6)}_B}{r_{scB}^{8}}\right)r^2 - \left(\frac{168C^{(12)}_B}{r_{scB}^{13}} - \frac{48C^{(6)}_B}{r_{scB}^{7}}\right)r
+ \frac{91C^{(12)}_B}{r_{scB}^{12}} - \frac{28C^{(6)}_B}{r_{scB}^{6}} \\
& -  \left[\left(\frac{78C^{(12)}_A}{r_{scA}^{14}} - \frac{21C^{(6)}_A}{r_{scA}^{8}}\right)r^2 - \left(\frac{168C^{(12)}_A}{r_{scA}^{13}} - \frac{48C^{(6)}_A}{r_{scA}^{7}}\right)r
+ \frac{91C^{(12)}_A}{r_{scA}^{12}} - \frac{28C^{(6)}_A}{r_{scA}^{6}} \right]\\
& +  \frac{14(\lambda-1)}{\lambda}\left[\left(\frac{13C^{(12)}_A}{r_{scA}^{14}} - \frac{2C^{(6)}_A}{r_{scA}^{8}}\right)r^2
- \left(\frac{26C^{(12)}_A}{r_{scA}^{13}} - \frac{4C^{(6)}_A}{r_{scA}^{7}}\right)r
+ \frac{13C^{(12)}_A}{r_{scA}^{12}} - \frac{2C^{(6)}_A}{r_{scA}^{6}}\right] \\
& +  \frac{14\lambda}{1-\lambda}\left[\left(\frac{13C^{(12)}_B}{r_{scB}^{14}} - \frac{2C^{(6)}_B}{r_{scB}^{8}}\right)r^2
- \left(\frac{26C^{(12)}_B}{r_{scB}^{13}} - \frac{4C^{(6)}_B}{r_{scB}^{7}}\right)r
+ \frac{13C^{(12)}_B}{r_{scB}^{12}} - \frac{2C^{(6)}_B}{r_{scB}^{6}}\right] \end{aligned}
\end{aligned}
$$

$\partial H/ \partial \lambda$ for Lennard-Jones potential, when
$r \geq r_{scLJ}$ is calculated as a standard hard-core contribution
to $\partial H/ \partial \lambda$:
$\frac{\partial{H}}{\partial{\lambda}} = V_{LJ}^B(r) - V_{LJ}^A(r)$.

### $\partial H/ \partial \lambda$ for Coulomb interactions

Here we provide the explicit expressions of
$\partial H/ \partial \lambda$ for Coulomb potential, when
$r<r_{scQ}<r_{cutoffQ}$. For simplicity, in the expression below we
use the notation $r_{scQ_A}=r_{scA}$ and $r_{scQ_B}=r_{scB}$.

$$
\begin{aligned}
\begin{aligned}
\frac{\partial{H}}{\partial{\lambda}} &= V_Q^B(r) - V_Q^A(r) + (1-\lambda)\frac{\partial{V_Q^A(r)}}{\partial{\lambda}} + \lambda\frac{\partial{V_Q^B(r)}}{\partial{\lambda}} \\
& =  \frac{q_{i}^Bq_{j}^B}{r_{scB}^3}r^2 - \frac{3q_i^Bq_j^B}{r_{scB}^2}r + \frac{3q^B_iq_j^B}{r_{scB}} \\
& -  \left[\frac{q_{i}^Aq_{j}^A}{r_{scA}^3}r^2 - \frac{3q_i^Aq_j^A}{r_{scA}^2}r + \frac{3q^A_iq_j^A}{r_{scA}}\right] \\
& +  \frac{\lambda-1}{2\lambda}\left[\frac{q_i^Aq_j^A}{r_{scA}^3}r^2 - \frac{2q_i^Aq_j^A}{r_{scA}^2}r + \frac{q_i^Aq_j^A}{r_{scA}}\right] \\
& +  \frac{\lambda}{2(1-\lambda)}\left[\frac{q_i^Bq_j^B}{r_{scB}^3}r^2 - \frac{2q_i^Bq_j^B}{r_{scB}^2}r + \frac{q_i^Bq_j^B}{r_{scB}}\right] \end{aligned}
\end{aligned}
$$

$\partial H/ \partial \lambda$ for Coulomb potential, when
$r < r_{scQ} \geq r_{cutoffQ}$ is calculated using the same expression
above by setting $r_{scA}=r_{cutoffQ}$ and $r_{scB}=r_{cutoffQ}$.

$\partial H/ \partial \lambda$ for Coulomb potential, when
$r \geq r_{scQ} < r_{cutoffQ}$ is calculated as a standard hard-core
contribution to $\partial H/ \partial \lambda$:
$\frac{\partial{H}}{\partial{\lambda}} = V_{Q}^B(r) - V_{Q}^A(r)$.
