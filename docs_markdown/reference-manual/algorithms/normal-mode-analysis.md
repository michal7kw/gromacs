# Normal-Mode Analysis

Normal-mode analysis `54 <refLevitt83>``56 <refBBrooks83b>` can be
performed using GROMACS, by diagonalization of
the mass-weighted Hessian $H$:

$$
\begin{aligned}
\begin{aligned}
R^T M^{-1/2} H M^{-1/2} R   &=& \mbox{diag}(\lambda_1,\ldots,\lambda_{3N})
\\
\lambda_i &=& (2 \pi \omega_i)^2\end{aligned}
\end{aligned}
$$

where $M$ contains the atomic masses, $R$ is a matrix that contains
the eigenvectors as columns, $\lambda_i$ are the eigenvalues and
$\omega_i$ are the corresponding frequencies.

First the Hessian matrix, which is a $3N \times 3N$ matrix where $N$
is the number of atoms, needs to be calculated:

$$
\begin{aligned}
H_{ij}  &=&     \frac{\partial^2 V}{\partial x_i \partial x_j}\end{aligned}
$$

where $x_i$ and $x_j$ denote the atomic x, y or z coordinates. In
practice, this equation is not used, but the Hessian is calculated
numerically from the force as:

$$
\begin{aligned}
\begin{aligned}
H_{ij} &=& -
\frac{f_i({\bf x}+h{\bf e}_j) - f_i({\bf x}-h{\bf e}_j)}{2h}
\\
f_i     &=& - \frac{\partial V}{\partial x_i}\end{aligned}
\end{aligned}
$$

where ${\bf e}_j$ is the unit vector in direction $j$. It should be
noted that for a usual normal-mode calculation, it is necessary to
completely minimize the energy prior to computation of the Hessian. The
tolerance required depends on the type of system, but a rough indication
is 0.001 kJ mol$^{-1}$. Minimization should be done with conjugate
gradients or L-BFGS in double precision.

A number of GROMACS programs are involved in
these calculations. First, the energy should be minimized using
`mdrun `. Then, `mdrun ` computes the Hessian.
**Note** that for generating the run input file, one should use the
minimized conformation from the full precision trajectory file, as the
structure file is not accurate enough. `gmx nmeig` does the
diagonalization and the sorting of the normal modes according to their
frequencies. Both `mdrun ` and `gmx nmeig` should be run in
double precision. The normal modes can be analyzed with the program
`gmx anaeig`. Ensembles of structures at any temperature and for any
subset of normal modes can be generated with `gmx nmens`. An overview of
normal-mode analysis and the related principal component analysis (see
sec. `covanal`) can be found in `57 <refHayward95b>`.
