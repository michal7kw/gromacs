# Averages and fluctuations

## Formulae for averaging

**Note:** this section was taken from ref `179 <refGunsteren94a>`.

When analyzing a MD trajectory averages $\left<x\right>$ and
fluctuations

$$
\left<(\Delta x)^2\right>^{{\frac{1}{2}}} ~=~ \left<[x-\left<x\right>]^2\right>^{{\frac{1}{2}}}
$$

of a quantity $x$ are to be computed. The variance $\sigma_x$ of a
series of N$_x$ values, {$x_i$}, can be computed from

$$
\sigma_x~=~ \sum_{i=1}^{N_x} x_i^2 ~-~  \frac{1}{N_x}\left(\sum_{i=1}^{N_x}x_i\right)^2
$$

Unfortunately this formula is numerically not very accurate, especially
when $\sigma_x^{{\frac{1}{2}}}$ is small compared to the values of
$x_i$. The following (equivalent) expression is numerically more
accurate

$$
\sigma_x ~=~ \sum_{i=1}^{N_x} [x_i  - \left<x\right>]^2
$$

with

$$
\left<x\right> ~=~ \frac{1}{N_x} \sum_{i=1}^{N_x} x_i
$$

Using `eqns. %s <eqnvar1>` and `%s <eqnvar2>` one has to go through the
series of $x_i$ values twice, once to determine $\left<x\right>$ and
again to compute $\sigma_x$, whereas `eqn. %s <eqnvar0>` requires only
one sequential scan of the series {$x_i$}. However, one may cast
`eqn. %s <eqnvar1>` in another form, containing partial sums, which
allows for a sequential update algorithm. Define the partial sum

$$
X_{n,m} ~=~ \sum_{i=n}^{m} x_i
$$

and the partial variance

$$
\sigma_{n,m} ~=~ \sum_{i=n}^{m}  \left[x_i - \frac{X_{n,m}}{m-n+1}\right]^2
$$

It can be shown that

$$
X_{n,m+k} ~=~  X_{n,m} + X_{m+1,m+k}
$$

and

$$
\begin{aligned}
\begin{aligned}
\sigma_{n,m+k} &=& \sigma_{n,m} + \sigma_{m+1,m+k} + \left[~\frac {X_{n,m}}{m-n+1} - \frac{X_{n,m+k}}{m+k-n+1}~\right]^2~* \nonumber\\
&& ~\frac{(m-n+1)(m+k-n+1)}{k}
\end{aligned}
\end{aligned}
$$

For $n=1$ one finds

$$
\sigma_{1,m+k} ~=~ \sigma_{1,m} + \sigma_{m+1,m+k}~+~
\left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^2~ \frac{m(m+k)}{k}
$$

and for $n=1$ and $k=1$ `eqn. %s <eqnvarpartial>` becomes

$$
\begin{aligned}
\begin{aligned}
\sigma_{1,m+1}  &=& \sigma_{1,m} +
\left[\frac{X_{1,m}}{m} - \frac{X_{1,m+1}}{m+1}\right]^2 m(m+1)\\
&=& \sigma_{1,m} +
\frac {[~X_{1,m} - m x_{m+1}~]^2}{m(m+1)}
\end{aligned}
\end{aligned}
$$

where we have used the relation

$$
X_{1,m+1} ~=~  X_{1,m} + x_{m+1}
$$

Using formulae `eqn. %s <eqnsimplevar0>` and `eqn. %s <eqnsimplevar1>`
the average

$$
\left<x\right> ~=~ \frac{X_{1,N_x}}{N_x}
$$

and the fluctuation

$$
\left<(\Delta x)^2\right>^{{\frac{1}{2}}} = \left[\frac {\sigma_{1,N_x}}{N_x}\right]^{{\frac{1}{2}}}
$$

can be obtained by one sweep through the data.

## Implementation

In GROMACS the instantaneous energies $E(m)$
are stored in the `energy file <edr>`, along with the values of
$\sigma_{1,m}$ and $X_{1,m}$. Although the steps are counted from 0,
for the energy and fluctuations steps are counted from 1. This means
that the equations presented here are the ones that are implemented. We
give somewhat lengthy derivations in this section to simplify checking
of code and equations later on.

### Part of a Simulation

It is not uncommon to perform a simulation where the first part, *e.g.*
100 ps, is taken as equilibration. However, the averages and
fluctuations as printed in the `log file <log>` are computed over the
whole simulation. The equilibration time, which is now part of the
simulation, may in such a case invalidate the averages and fluctuations,
because these numbers are now dominated by the initial drift towards
equilibrium.

Using `eqns. %s <eqnXpartial>` and `%s <eqnvarpartial>` the average and
standard deviation over part of the trajectory can be computed as:

$$
\begin{aligned}
\begin{aligned}
X_{m+1,m+k}     &=& X_{1,m+k} - X_{1,m}                 \\
\sigma_{m+1,m+k} &=& \sigma_{1,m+k}-\sigma_{1,m} - \left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^{2}~ \frac{m(m+k)}{k}\end{aligned}
\end{aligned}
$$

or, more generally (with $p \geq 1$ and $q \geq p$):

$$
\begin{aligned}
\begin{aligned}
X_{p,q}         &=&     X_{1,q} - X_{1,p-1}     \\
\sigma_{p,q}    &=&     \sigma_{1,q}-\sigma_{1,p-1} - \left[~\frac{X_{1,p-1}}{p-1} - \frac{X_{1,q}}{q}~\right]^{2}~ \frac{(p-1)q}{q-p+1}\end{aligned}
\end{aligned}
$$

**Note** that implementation of this is not entirely trivial, since
energies are not stored every time step of the simulation. We therefore
have to construct $X_{1,p-1}$ and $\sigma_{1,p-1}$ from the
information at time $p$ using `eqns. %s <eqnsimplevar0>` and
`%s <eqnsimplevar1>`:

$$
\begin{aligned}
\begin{aligned}
X_{1,p-1}       &=&     X_{1,p} - x_p   \\
\sigma_{1,p-1}  &=&     \sigma_{1,p} -  \frac {[~X_{1,p-1} - (p-1) x_{p}~]^2}{(p-1)p}\end{aligned}
\end{aligned}
$$

### Combining two simulations

Another frequently occurring problem is, that the fluctuations of two
simulations must be combined. Consider the following example: we have
two simulations (A) of $n$ and (B) of $m$ steps, in which the second
simulation is a continuation of the first. However, the second
simulation starts numbering from 1 instead of from $n+1$. For the
partial sum this is no problem, we have to add $X_{1,n}^A$ from run A:

$$
X_{1,n+m}^{AB} ~=~ X_{1,n}^A + X_{1,m}^B
$$

When we want to compute the partial variance from the two components we
have to make a correction $\Delta\sigma$:

$$
\sigma_{1,n+m}^{AB} ~=~ \sigma_{1,n}^A + \sigma_{1,m}^B +\Delta\sigma
$$

if we define $x_i^{AB}$ as the combined and renumbered set of data
points we can write:

$$
\sigma_{1,n+m}^{AB} ~=~ \sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2
$$

and thus

$$
\sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2  ~=~
\sum_{i=1}^{n}  \left[x_i^{A} - \frac{X_{1,n}^{A}}{n}\right]^2  +
\sum_{i=1}^{m}  \left[x_i^{B} - \frac{X_{1,m}^{B}}{m}\right]^2  +\Delta\sigma
$$

or

$$
\begin{aligned}
\begin{aligned}
\sum_{i=1}^{n+m}  \left[(x_i^{AB})^2 - 2 x_i^{AB}\frac{X^{AB}_{1,n+m}}{n+m} + \left(\frac{X^{AB}_{1,n+m}}{n+m}\right)^2  \right] &-& \nonumber \\
\sum_{i=1}^{n}  \left[(x_i^{A})^2 - 2 x_i^{A}\frac{X^A_{1,n}}{n} + \left(\frac{X^A_{1,n}}{n}\right)^2  \right] &-& \nonumber \\
\sum_{i=1}^{m}  \left[(x_i^{B})^2 - 2 x_i^{B}\frac{X^B_{1,m}}{m} + \left(\frac{X^B_{1,m}}{m}\right)^2  \right] &=& \Delta\sigma\end{aligned}
\end{aligned}
$$

all the $x_i^2$ terms drop out, and the terms independent of the
summation counter $i$ can be simplified:

$$
\begin{aligned}
\begin{aligned}
\frac{\left(X^{AB}_{1,n+m}\right)^2}{n+m} \,-\,
\frac{\left(X^A_{1,n}\right)^2}{n} \,-\,
\frac{\left(X^B_{1,m}\right)^2}{m} &-& \nonumber \\
2\,\frac{X^{AB}_{1,n+m}}{n+m}\sum_{i=1}^{n+m}x_i^{AB} \,+\,
2\,\frac{X^{A}_{1,n}}{n}\sum_{i=1}^{n}x_i^{A} \,+\,
2\,\frac{X^{B}_{1,m}}{m}\sum_{i=1}^{m}x_i^{B} &=& \Delta\sigma\end{aligned}
\end{aligned}
$$

we recognize the three partial sums on the second line and use
`eqn. %s <eqnpscomb>` to obtain:

$$
\Delta\sigma ~=~ \frac{\left(mX^A_{1,n} - nX^B_{1,m}\right)^2}{nm(n+m)}
$$

if we check this by inserting $m=1$ we get back
`eqn. %s <eqnsimplevar0>`

### Summing energy terms

The `gmx energy ` program can also sum energy terms into
one, *e.g.* potential + kinetic = total. For the partial averages this
is again easy if we have $S$ energy components $s$:

$$
X_{m,n}^S ~=~ \sum_{i=m}^n \sum_{s=1}^S x_i^s ~=~ \sum_{s=1}^S \sum_{i=m}^n x_i^s ~=~ \sum_{s=1}^S X_{m,n}^s
$$

For the fluctuations it is less trivial again, considering for example
that the fluctuation in potential and kinetic energy should cancel.
Nevertheless we can try the same approach as before by writing:

$$
\sigma_{m,n}^S ~=~ \sum_{s=1}^S \sigma_{m,n}^s + \Delta\sigma
$$

if we fill in `eqn. %s <eqnsigma>`:

$$
\sum_{i=m}^n \left[\left(\sum_{s=1}^S x_i^s\right) - \frac{X_{m,n}^S}{m-n+1}\right]^2 ~=~
\sum_{s=1}^S \sum_{i=m}^n \left[\left(x_i^s\right) - \frac{X_{m,n}^s}{m-n+1}\right]^2 + \Delta\sigma
$$

which we can expand to:

$$
\begin{aligned}
\begin{aligned}
&~&\sum_{i=m}^n \left[\sum_{s=1}^S (x_i^s)^2 + \left(\frac{X_{m,n}^S}{m-n+1}\right)^2 -2\left(\frac{X_{m,n}^S}{m-n+1}\sum_{s=1}^S x_i^s + \sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'} \right)\right]    \nonumber \\
&-&\sum_{s=1}^S \sum_{i=m}^n \left[(x_i^s)^2 - 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}
\end{aligned}
$$

the terms with $(x_i^s)^2$ cancel, so that we can simplify to:

$$
\begin{aligned}
\begin{aligned}
&~&\frac{\left(X_{m,n}^S\right)^2}{m-n+1} -2 \frac{X_{m,n}^S}{m-n+1}\sum_{i=m}^n\sum_{s=1}^S x_i^s -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, -        \nonumber \\
&~&\sum_{s=1}^S \sum_{i=m}^n \left[- 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}
\end{aligned}
$$

or

$$
-\frac{\left(X_{m,n}^S\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +  \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma
$$

If we now expand the first term using `eqn. %s <eqnsumterms>` we obtain:

$$
-\frac{\left(\sum_{s=1}^SX_{m,n}^s\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +      \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma
$$

which we can reformulate to:

$$
-2\left[\sum_{s=1}^S \sum_{s'=s+1}^S X_{m,n}^s X_{m,n}^{s'}\,+\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\right] ~=~\Delta\sigma
$$

or

$$
-2\left[\sum_{s=1}^S X_{m,n}^s \sum_{s'=s+1}^S X_{m,n}^{s'}\,+\,\sum_{s=1}^S \sum_{i=m}^nx_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma
$$

which gives

$$
-2\sum_{s=1}^S \left[X_{m,n}^s \sum_{s'=s+1}^S \sum_{i=m}^n x_i^{s'}\,+\,\sum_{i=m}^n x_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma
$$

Since we need all data points $i$ to evaluate this, in general this is
not possible. We can then make an estimate of $\sigma_{m,n}^S$ using
only the data points that are available using the left hand side of
`eqn. %s <eqnsigmaterms>`. While the average can be computed using all
time steps in the simulation, the accuracy of the fluctuations is thus
limited by the frequency with which energies are saved. Since this can
be easily done with a program such as `xmgr` this is not built-in in
GROMACS.
