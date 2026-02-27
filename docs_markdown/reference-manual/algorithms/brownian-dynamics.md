# Brownian Dynamics

In the limit of high friction, stochastic dynamics reduces to Brownian
dynamics, also called position Langevin dynamics. This applies to
over-damped systems, *i.e.* systems in which the inertia effects are
negligible. The equation is

$$
{{\mbox{d}}\mathbf{r}_i \over {\mbox{d}}t} = \frac{1}{\gamma_i} \mathbf{F}_i(\mathbf{r}) + {\stackrel{\circ}{\mathbf{r}}}_i
$$

where $\gamma_i$ is the friction coefficient $[\mbox{amu/ps}]$ and
${\stackrel{\circ}{\mathbf{r}}}_i(t)$ is a noise process with
$\langle {\stackrel{\circ}{r}}_i(t) {\stackrel{\circ}{r}}_j(t+s) \rangle = 2 \delta(s) \delta_{ij} k_B T / \gamma_i$.
In GROMACS the equations are integrated with a
simple, explicit scheme

$$
\mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) +
{\Delta t \over \gamma_i} \mathbf{F}_i(\mathbf{r}(t))
+ \sqrt{2 k_B T {\Delta t \over \gamma_i}}\, {\mathbf{r}^G}_i,
$$

where ${\mathbf{r}^G}_i$ is Gaussian distributed noise with
$\mu = 0$, $\sigma = 1$. The friction coefficients $\gamma_i$ can
be chosen the same for all particles or as $\gamma_i = m_i\,\gamma_i$,
where the friction constants $\gamma_i$ can be different for different
groups of atoms. Because the system is assumed to be over-damped, large
timesteps can be used. LINCS should be used for the constraints since
SHAKE will not converge for large atomic displacements. BD can be
activated by using `integrator=bd` and the simulations are run using the
`mdrun ` program.

In BD there are no velocities, so there is also no kinetic energy. Still
`gmx mdrun` will report a kinetic energy and temperature based on atom
displacements per step $\Delta x$. This can be used to judge the
quality of the integration. A too high temperature is an indication that
the time step chosen is too large. The formula for the kinetic energy
term reported is:

$$
\frac{1}{2} \sum_i \frac{\gamma_i \Delta x_i^2}{2 \, \Delta t}
$$
