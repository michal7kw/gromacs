# Shell molecular dynamics

GROMACS can simulate polarizability using the
shell model of Dick and Overhauser `43 <refDick58>`. In such models a
shell particle representing the electronic degrees of freedom is
attached to a nucleus by a spring. The potential energy is minimized
with respect to the shell position at every step of the simulation (see
below). Successful applications of shell models in
GROMACS have been published for $N_2$
`44 <refJordan95>` and water `45 <refMaaren2001a>`.

## Optimization of the shell positions

The force $\mathbf{F}_S$ on a shell particle $S$ can be decomposed
into two components

$$
\mathbf{F}_S ~=~ \mathbf{F}_{bond} + \mathbf{F}_{nb}
$$

where $\mathbf{F}_{bond}$ denotes the component representing the
polarization energy, usually represented by a harmonic potential and
$\mathbf{F}_{nb}$ is the sum of Coulomb and van der Waals
interactions. If we assume that $\mathbf{F}_{nb}$ is almost constant
we can analytically derive the optimal position of the shell, i.e. where
$\mathbf{F}_S = 0$. If we have the shell S connected to atom A we have

$$
\mathbf{F}_{bond} ~=~ k_b \left( \mathbf{x}_S - \mathbf{x}_A\right).
$$

In an iterative solver, we have positions $\mathbf{x}_S(n)$ where
$n$ is the iteration count. We now have at iteration $n$

$$
\mathbf{F}_{nb} ~=~ \mathbf{F}_S - k_b \left( \mathbf{x}_S(n) - \mathbf{x}_A\right)
$$

and the optimal position for the shells $x_S(n+1)$ thus follows from

$$
\mathbf{F}_S - k_b \left( \mathbf{x}_S(n) - \mathbf{x}_A\right) + k_b \left( \mathbf{x}_S(n+1) - \mathbf{x}_A\right) = 0
$$

if we write

$$
\Delta \mathbf{x}_S = \mathbf{x}_S(n+1) - \mathbf{x}_S(n)
$$

we finally obtain

$$
\Delta \mathbf{x}_S = \mathbf{F}_S/k_b
$$

which then yields the algorithm to compute the next trial in the
optimization of shell positions

$$
\mathbf{x}_S(n+1) ~=~ \mathbf{x}_S(n) + \mathbf{F}_S/k_b.
$$
