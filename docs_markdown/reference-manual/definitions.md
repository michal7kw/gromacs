# Definitions and Units

## Notation

The following conventions for mathematical typesetting are used
throughout this document:

| Item          | Notation    | Example            |
|---------------|-------------|--------------------|
| Vector        | Bold italic | ${\mathbf{r}_i}$ |
| Vector Length | Italic      | $r_i$            |

We define the *lowercase* subscripts $i$, $j$, $k$ and $l$ to
denote particles: $\mathbf{r}_i$ is the *position vector* of particle
$i$, and using this notation:

$$
\begin{aligned}
\begin{aligned}
\mathbf{r}_{ij} = \mathbf{r}_j-\mathbf{r}_i \\
r_{ij}          = | \mathbf{r}_{ij} |       \end{aligned}
\end{aligned}
$$

The force on particle $i$ is denoted by $\mathbf{F}_i$ and

$$
\mathbf{F}_{ij} = \mbox{force on $i$ exerted by $j$}
$$

## MD units

GROMACS uses a consistent set of units that
produce values in the vicinity of unity for most relevant molecular
quantities. Let us call them *MD units*. The basic units in this system
are nm, ps, K, electron charge (e) and atomic mass unit (u), see
`Table %s ` The values used in
GROMACS are taken from the CODATA
Internationally recommended 2010 values of fundamental physical
constants (see [NIST homepage](http://nist.gov)).

**Basic units used in |Gromacs|**

| Quantity | Symbol | Unit |
| --- | --- | --- |
| length | r | $nm=10−9 m$ |
| mass | m | u (unified atomic mass unit) = $1.660 538 921 × 10−27 kg$ |
| time | t | $ps=10−12 s$ |
| charge | q | e = elementary charge = $1.602 176 565 × 10−19 C$ |
| temperature | T | K |

Consistent with these units are a set of derived units, given in
`Table %s `

| Quantity | Symbol | Unit |
|----|----|----|
| energy | $E,V$ | $\mathrm{kJ~mol}^{-1}$ |
| Force | $\mathbf{F}$ | $\mathrm{kJ~mol}^{-1}~\mathrm{nm}^{-1}$ |
| pressure | $p$ | bar |
| velocity | $v$ | $\mathrm{nm~ps}^{-1} = 1000\mathrm{~m~s}^{-1}$ |
| dipole moment | $\mu$ | $\mathrm{e\ nm}$ |
| electric potential | $\Phi$ | $\mathrm{kJ~mol}^{-1}\mathrm{~e}^{-1} =$ $0.010\,364\,269\,19$ Volt |
| electric field | $E$ | $\mathrm{kJ~mol}^{-1}\mathrm{~nm}^{-1}\ \mathrm{e}^{-1} =$ $1.036\,426\,919 \times 10^7\mathrm{~V m}^{-1}$ |

Derived units. Note that an additional conversion factor of 10$^{28}$
a.m.u ($\approx$ 16.6) is applied to get bar instead of internal MD
units in the energy and log files

The **electric conversion factor**
$f=\frac{1}{4 \pi \varepsilon_o}={138.935\,458}$
$\mathrm{kJ}~\mathrm{mol}^{-1}\mathrm{nm}~\mathrm{ e}^{-2}$. It
relates the mechanical quantities to the electrical quantities as in

$$
V = f \frac{q^2}{r} \mbox{\ \ or\ \ } F = f \frac{q^2}{r^2}
$$

Electric potentials $\Phi$ and electric fields $\mathbf{E}$ are
intermediate quantities in the calculation of energies and forces. They
do not occur inside GROMACS. If they are used
in evaluations, there is a choice of equations and related units. We
strongly recommend following the usual practice of including the factor
$f$ in expressions that evaluate $\Phi$ and $\mathbf{E}$:

$$
\begin{aligned}
\begin{aligned}
\Phi(\mathbf{r}) = f \sum_j \frac{q_j}{| \mathbf{r}-\mathbf{r}_j | } \\
\mathbf{E}(\mathbf{r}) = f \sum_j q_j \frac{(\mathbf{r}-\mathbf{r}_j)}{| \mathbf{r}-\mathbf{r}_j| ^3}\end{aligned}
\end{aligned}
$$

With these definitions, $q\Phi$ is an energy and $q\mathbf{E}$ is a
force. The units are those given in `Table %s `
about 10 mV for potential. Thus, the potential of an electronic charge
at a distance of 1 nm equals $f \approx 140$ units $\approx 1.4$ V.
(exact value: $1.439\,964\,5$ V)

**Note** that these units are mutually consistent; changing any of the
units is likely to produce inconsistencies and is therefore *strongly
discouraged*! In particular: if Å are used instead of nm, the unit of
time changes to 0.1 ps. If $\mathrm{kcal}~\mathrm{mol}^{-1}$ (= 4.184
$\mathrm{kJ~mol}^{-1}$) is used instead of $\mathrm{kJ~mol}^{-1}$
for energy, the unit of time becomes 0.488882 ps and the unit of
temperature changes to 4.184 K. But in both cases all electrical
energies go wrong, because they will still be computed in
$\mathrm{kJ~mol}^{-1}$, expecting nm as the unit of length. Although
careful rescaling of charges may still yield consistency, it is clear
that such confusions must be rigidly avoided.

In terms of the MD units, the usual physical constants take on different
values (see `Table %s `). All quantities are per mol
rather than per molecule. There is no distinction between Boltzmann’s
constant $k$ and the gas constant $R$: their value is
$0.008\,314\,462\,1\mathrm{kJ~mol}^{-1} \mathrm{K}^{-1}$.

| Symbol | Name | Value |
|----|----|----|
| $N_{AV}$ | Avogadro's number | $6.022\,141\,29\times 10^{23}~\mathrm{mol}^{-1}$ |
| $R$ | gas constant | $8.314\,462\,1\times 10^{-3}~\mathrm{kJ~mol}^{-1}~\mathrm{K}^{-1}$ |
| $k_B$ | Boltzmann's constant | *idem* |
| $h$ | Planck's constant | $0.399\,031\,271~\mathrm{kJ~mol}^{-1}~\mathrm{ps}$ |
| $\hbar$ | Dirac's constant | $0.063\,507\,799\,3~\mathrm{kJ~mol}^{-1}~\mathrm{ps}$ |
| $c$ | velocity of light | $299\,792.458~\mathrm{nm~ps}^{-1}$ |

Some Physical Constants

## Reduced units

When simulating Lennard-Jones (LJ) systems, it might be advantageous to
use reduced units (*i.e.*, setting
$\epsilon_{ii}=\sigma_{ii}=m_i=k_B=1$ for one type of atoms). This is
possible. When specifying the input in reduced units, the output will
also be in reduced units. The one exception is the *temperature*, which
is expressed in $0.008\,314\,462\,1$ reduced units. This is a
consequence of using Boltzmann’s constant in the evaluation of
temperature in the code. Thus not $T$, but $k_BT$, is the reduced
temperature. A GROMACS temperature $T=1$
means a reduced temperature of $0.008\ldots$ units; if a reduced
temperature of 1 is required, the GROMACS
temperature should be $120.272\,36$.

In `Table %s ` quantities are given for LJ potentials:

$$
V_{LJ} = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]
$$

| Quantity    | Symbol     | Relation to SI                     |
|-------------|------------|------------------------------------|
| Length      | r$^*$    | r$\sigma^{-1}$                   |
| Mass        | m$^*$    | m M$^{-1}$                       |
| Time        | t$^*$    | t$\sigma^{-1}~\sqrt{\epsilon/M}$ |
| Temperature | T$^*$    | k$_B\mathrm{T}~\epsilon^{-1}$    |
| Energy      | E$^*$    | E$\epsilon^{-1}$                 |
| Force       | F$^*$    | F$\sigma~\epsilon^{-1}$          |
| Pressure    | P$^*$    | P$\sigma ^3 \epsilon^{-1}$       |
| Velocity    | v$^*$    | v$\sqrt{M/\epsilon}$             |
| Density     | $\rho^*$ | N$\sigma ^3~V^{-1}$              |

Reduced Lennard-Jones quantities

## Mixed or Double precision

GROMACS can be compiled in either mixed or
double precision. Documentation of previous
GROMACS versions referred to *single
precision*, but the implementation has made selective use of double
precision for many years. Using single precision for all variables would
lead to a significant reduction in accuracy. Although in *mixed
precision* all state vectors, i.e. particle coordinates, velocities and
forces, are stored in single precision, critical variables are double
precision. A typical example of the latter is the virial, which is a sum
over all forces in the system, which have varying signs. In addition, in
many parts of the code we managed to avoid double precision for
arithmetic, by paying attention to summation order or reorganization of
mathematical expressions. The default configuration uses mixed
precision, but it is easy to turn on double precision by adding the
option `-DGMX_DOUBLE=on` to `cmake`. Double precision will be 20 to 100%
slower than mixed precision depending on the architecture you are
running on. Double precision will use somewhat more memory and run
input, energy and full-precision trajectory files will be almost twice
as large.

The energies in mixed precision are accurate up to the last decimal, the
last one or two decimals of the forces are non-significant. The virial
is less accurate than the forces, since the virial is only one order of
magnitude larger than the size of each element in the sum over all atoms
(sec. `virial`). In most cases this is not really a problem, since the
fluctuations in the virial can be two orders of magnitude larger than
the average. Using cut-offs for the Coulomb interactions cause large
errors in the energies, forces, and virial. Even when using a
reaction-field or lattice sum method, the errors are larger than, or
comparable to, the errors due to the partial use of single precision.
Since MD is chaotic, trajectories with very similar starting conditions
will diverge rapidly, the divergence is faster in mixed precision than
in double precision.

For most simulations, mixed precision is accurate enough. In some cases
double precision is required to get reasonable results:

- normal mode analysis, for the conjugate gradient or l-bfgs
  minimization and the calculation and diagonalization of the Hessian
- long-term energy conservation, especially for large systems
