# Energy Minimization

Energy minimization in GROMACS can be done
using steepest descent, conjugate gradients, or l-bfgs (limited-memory
Broyden-Fletcher-Goldfarb-Shanno quasi-Newtonian minimizer...we prefer
the abbreviation). Whether to use EM, and which algorithm to use, is
specified via the `integrator` setting of the `mdrun `
program.

## Steepest Descent

Although steepest descent is certainly not the most efficient algorithm
for searching, it is robust and easy to implement.

We define the vector $\mathbf{r}$ as the vector of all $3N$
coordinates. Initially a maximum displacement $h_0$ (*e.g.* 0.01 nm)
must be given.

First the forces $\mathbf{F}$ and potential energy are calculated. New
positions are calculated by

$$
\mathbf{r}_{n+1} =  \mathbf{r}_n + \frac{\mathbf{F}_n}{\max (|\mathbf{F}_n|)} h_n,

$$

where $h_n$ is the maximum displacement and $\mathbf{F}_n$ is the
force, or the negative gradient of the potential $V$. The notation
$\max
(|\mathbf{F}_n|)$ means the largest scalar force on any atom. The
forces and energy are again computed for the new positions

If ($V_{n+1} < V_n$) the new positions are accepted and
$h_{n+1} = 1.2  h_n$.  
If ($V_{n+1} \geq V_n$) the new positions are rejected and
$h_n = 0.2 h_n$.

The algorithm stops when either a user-specified number of force
evaluations has been performed (*e.g.* 100), or when the maximum of the
absolute values of the force (gradient) components is smaller than a
specified value $\epsilon$. Since force truncation produces some noise
in the energy evaluation, the stopping criterion should not be made too
tight to avoid endless iterations. A reasonable value for $\epsilon$
can be estimated from the root mean square force $f$ a harmonic
oscillator would exhibit at a temperature $T$. This value is

$$
f = 2 \pi \nu \sqrt{ 2mkT},
$$

where $\nu$ is the oscillator frequency, $m$ the (reduced) mass, and
$k$ Boltzmann’s constant. For a weak oscillator with a wave number of
100 cm$^{-1}$ and a mass of 10 atomic units, at a temperature of 1 K,
$f=7.7$ kJ mol$^{-1}$ nm$^{-1}$. A value for $\epsilon$ between
1 and 10 is acceptable.

## Conjugate Gradient

Conjugate gradient is slower than steepest descent in the early stages
of the minimization, but becomes more efficient closer to the energy
minimum. The parameters and stop criterion are the same as for steepest
descent. The most common use case for conjugate gradient is minimization
prior to a normal-mode analysis, which requires very small forces. For
most other purposes steepest descent is efficient enough.

## L-BFGS

The original BFGS algorithm works by successively creating better
approximations of the inverse Hessian matrix, and moving the system to
the currently estimated minimum. The memory requirements for this are
proportional to the square of the number of particles, so it is not
practical for large systems like biomolecules. Instead, we use the
L-BFGS algorithm of Nocedal  `52 <refByrd95a>`, `53 <refZhu97a>`, which
approximates the inverse Hessian by a fixed number of corrections from
previous steps. This sliding-window technique is almost as efficient as
the original method, but the memory requirements are much lower -
proportional to the number of particles multiplied with the correction
steps. In practice we have found it to converge faster than conjugate
gradients, but due to the correction steps it is not yet parallelized.
It is also noteworthy that switched or shifted interactions usually
improve the convergence, since sharp cut-offs mean the potential
function at the current coordinates is slightly different from the
previous steps used to build the inverse Hessian approximation.
