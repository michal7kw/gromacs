# Some implementation details

In this chapter we will present some implementation details. This is far
from complete, but we deemed it necessary to clarify some things that
would otherwise be hard to understand.

## Single Sum Virial in GROMACS

The virial $\Xi$ can be written in full tensor form as:

$$
\Xi~=~-\frac{1}{2}~\sum_{i < j}^N~\mathbf{r}_{ij}\otimes\mathbf{F}_{ij}
$$

where $\otimes$ denotes the *direct product* of two vectors.[^1] When
this is computed in the inner loop of an MD program 9 multiplications
and 9 additions are needed.[^2]

Here it is shown how it is possible to extract the virial calculation
from the inner loop `177 <refBekker93b>`.

### Virial

In a system with periodic boundary conditions, the periodicity must be
taken into account for the virial:

$$
\Xi~=~-\frac{1}{2}~\sum_{i < j}^{N}~\mathbf{r}_{ij}^n\otimes\mathbf{F}_{ij}
$$

where $\mathbf{r}_{ij}^n$ denotes the distance vector of the *nearest
image* of atom $i$ from atom $j$. In this definition we add a *shift
vector* $\delta_i$ to the position vector $\mathbf{r}_i$ of atom
$i$. The difference vector $\mathbf{r}_{ij}^n$ is thus equal to:

$$
\mathbf{r}_{ij}^n~=~\mathbf{r}_i+\delta_i-\mathbf{r}_j
$$

or in shorthand:

$$
\mathbf{r}_{ij}^n~=~\mathbf{r}_i^n-\mathbf{r}_j
$$

In a triclinic system, there are 27 possible images of $i$; when a
truncated octahedron is used, there are 15 possible images.

### Virial from non-bonded forces

Here the derivation for the single sum virial in the *non-bonded force*
routine is given. There are a couple of considerations that are special
to GROMACS that we take into account:

- When calculating short-range interactions, we apply the *minimum image
  convention* and only consider the closest image of each neighbor - and
  in particular we never allow interactions between a particle and any
  of its periodic images. For all the equations below, this means
  $i \neq j$.
- In general, either the $i$ or $j$ particle might be shifted to a
  neighbor cell to get the closest interaction (shift $\delta_{ij}$).
  However, with minimum image convention there can be at most 27
  different shifts for particles in the central cell, and for typical
  (very short-ranged) biomolecular interactions there are typically only
  a few different shifts involved for each particle, not to mention that
  each interaction can only be present for one shift.
- For the GROMACS nonbonded interactions we
  use this to split the neighborlist of each $i$ particle into
  multiple separate lists, where each list has a constant shift
  $\delta_i$ for the $i$ partlcle. We can represent this as a sum
  over shifts (for which we use index $s$), with the constraint that
  each particle interaction can only contribute to one of the terms in
  this sum, and the shift is no longer dependent on the $j$ particles.
  For any sum that does not contain complex dependence on $s$, this
  means the sum trivially reduces to just the sum over $i$ and/or
  $j$.
- To simplify some of the sums, we replace sums over $j<i$ with double
  sums over all particles (remember, $i \neq j$) and divide by 2.

Starting from the above definition of the virial, we then get

$$
\begin{aligned}
\begin{aligned}
\Xi
&~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~{\mathbf r}^n_{ij} \otimes {\mathbf F}_{ij} \nonumber \\
&~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij,s} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{i=}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( \left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} -{\mathbf r}_j \otimes {\mathbf F}_{ij,s} \right) \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij,s} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ij} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} - {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ji} \nonumber \\
&~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{i,s} - {\frac{1}{4}}\sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j}  \nonumber \\
&~=~&-{\frac{1}{4}}~\left(\sum_{i=1}^{N}~{\mathbf r}_i  \otimes {\mathbf F}_{i} + \sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j} \right) - {\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N} \delta_{i,s} \otimes {\mathbf F}_{i,s}  \nonumber \\
&~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N}~\delta_{i,s} \otimes {\mathbf F}_{i,s} \nonumber \\
&~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\delta_{s} \otimes {\mathbf F}_{s} \nonumber \\
&~=~&\Xi_0 + \Xi_1\end{aligned}
\end{aligned}
$$

In the second-last stage, we have used the property that each shift
vector itself does not depend on the coordinates of particle $i$, so
it is possible to sum up all forces corresponding to each shift vector
(in the nonbonded kernels), and then just use a sum over the different
shift vectors outside the kernels. We have also used

$$
\begin{aligned}
\begin{aligned}
\mathbf{F}_i&~=~&\sum_{j=1}^N~\mathbf{F}_{ij} \\
\mathbf{F}_j&~=~&\sum_{i=1}^N~\mathbf{F}_{ji} \end{aligned}
\end{aligned}
$$

which is the total force on $i$ with respect to $j$. Because we use
Newton’s Third Law:

$$
\mathbf{F}_{ij}~=~-\mathbf{F}_{ji}
$$

we must, in the implementation, double the term containing the shift
$\delta_i$. Similarly, in a few places we have summed the
shift-dependent force over all shifts to come up with the total force
per interaction or particle.

This separates the total virial $\Xi$ into a component $\Xi_0$ that
is a single sum over particles, and a second component $\Xi_1$ that
describes the influence of the particle shifts, and that is only a sum
over the different shift vectors.

### The intra-molecular shift (mol-shift)

For the bonded forces and SHAKE it is possible to make a *mol-shift*
list, in which the periodicity is stored. We simple have an array mshift
in which for each atom an index in the shiftvec array is stored.

The algorithm to generate such a list can be derived from graph theory,
considering each particle in a molecule as a bead in a graph, the bonds
as edges.

1.  Represent the bonds and atoms as bidirectional graph
2.  Make all atoms white
3.  Make one of the white atoms black (atom $i$) and put it in the
    central box
4.  Make all of the neighbors of $i$ that are currently white, gray
5.  Pick one of the gray atoms (atom $j$), give it the correct
    periodicity with respect to any of its black neighbors and make it
    black
6.  Make all of the neighbors of $j$ that are currently white, gray
7.  If any gray atom remains, go to \[5\]
8.  If any white atom remains, go to \[3\]

Using this algorithm we can

- optimize the bonded force calculation as well as SHAKE
- calculate the virial from the bonded forces in the single sum method
  again

Find a representation of the bonds as a bidirectional graph.

### Virial from Covalent Bonds

Since the covalent bond force gives a contribution to the virial, we
have:

$$
\begin{aligned}
\begin{aligned}
b            &~=~& \|\mathbf{r}_{ij}^n\|                  \\
V_b          &~=~& \frac{1}{2} k_b(b-b_0)^2               \\
\mathbf{F}_i &~=~& -\nabla V_b                            \\
&~=~& k_b(b-b_0)\frac{\mathbf{r}_{ij}^n}{b}  \\
\mathbf{F}_j &~=~& -\mathbf{F}_i\end{aligned}
\end{aligned}
$$

The virial contribution from the bonds then is:

$$
\begin{aligned}
\begin{aligned}
\Xi_b &~=~& -\frac{1}{2}(\mathbf{r}_i^n\otimes\mathbf{F}_i~+~\mathbf{r}_j\otimes\mathbf{F}_j) \\
&~=~& -\frac{1}{2}\mathbf{r}_{ij}^n\otimes\mathbf{F}_i\end{aligned}
\end{aligned}
$$

### Virial from SHAKE

An important contribution to the virial comes from shake. Satisfying the
constraints a force **G** that is exerted on the particles “shaken.” If
this force does not come out of the algorithm (as in standard SHAKE) it
can be calculated afterward (when using *leap-frog*) by:

$$
\begin{aligned}
\begin{aligned}
\Delta\mathbf{r}_i&~=~&{\mathbf{r}_i}(t+{\Delta t})-
[\mathbf{r}_i(t)+{\bf v}_i(t-\frac{{\Delta t}}{2}){\Delta t}+\frac{\mathbf{F}_i}{m_i}{\Delta t}^2]    \\
{\bf G}_i&~=~&\frac{m_i{\Delta}{\mathbf{r}_i}}{{\Delta t}^2i}\end{aligned}
\end{aligned}
$$

This does not help us in the general case. Only when no periodicity is
needed (like in rigid water) this can be used, otherwise we must add the
virial calculation in the inner loop of SHAKE.

When it *is* applicable the virial can be calculated in the single sum
way:

$$
\Xi~=~-\frac{1}{2}\sum_i^{N_c}~\mathbf{r}_i\otimes\mathbf{F}_i
$$

where $N_c$ is the number of constrained atoms.

## Optimizations

Here we describe some of the algorithmic optimizations used in
GROMACS, apart from parallelism.

### Inner Loops for Water

GROMACS uses special inner loops to calculate
non-bonded interactions for water molecules with other atoms, and yet
another set of loops for interactions between pairs of water molecules.
There highly optimized loops for two types of water models. For three
site models similar to SPC `80 <refBerendsen81>`, *i.e.*:

1.  There are three atoms in the molecule.
2.  The whole molecule is a single charge group.
3.  The first atom has Lennard-Jones (sec. `lj`) and Coulomb
    (sec. `coul`) interactions.
4.  Atoms two and three have only Coulomb interactions, and equal
    charges.

These loops also works for the SPC/E `178 <refBerendsen87>` and
TIP3P `128 <refJorgensen83>` water models. And for four site water
models similar to TIP4P `128 <refJorgensen83>`:

1.  There are four atoms in the molecule.
2.  The whole molecule is a single charge group.
3.  The first atom has only Lennard-Jones (sec. `lj`) interactions.
4.  Atoms two and three have only Coulomb (sec. `coul`) interactions,
    and equal charges.
5.  Atom four has only Coulomb interactions.

The benefit of these implementations is that there are more
floating-point operations in a single loop, which implies that some
compilers can schedule the code better. However, it turns out that even
some of the most advanced compilers have problems with scheduling,
implying that manual tweaking is necessary to get optimum performance.
This may include common-sub-expression elimination, or moving code
around.

[^1]: Note that some derivations, an alternative notation
    $\xi_{\mathrm{alt}} = v_{\xi} = p_{\xi}/Q$ is used.

[^2]: The calculation of Lennard-Jones and Coulomb forces is about 50
    floating point operations.
