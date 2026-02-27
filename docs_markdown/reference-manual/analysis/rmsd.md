# Root mean square deviations in structure

`gmx rms `, `gmx rmsdist `  
The *root mean square deviation* ($RMSD$) of certain atoms in a
molecule with respect to a reference structure can be calculated with
the program `gmx rms ` by least-square fitting the structure to
the reference structure ($t_2 = 0$) and subsequently calculating the
$RMSD$ (`eqn. %s <eqnrmsd>`).

$$
RMSD(t_1,t_2) ~=~ \left[\frac{1}{M} \sum_{i=1}^N m_i \|{\bf r}_i(t_1)-{\bf r}_i(t_2)\|^2 \right]^{\frac{1}{2}}

$$

where $M = \sum_{i=1}^N m_i$ and ${\bf r}_i(t)$ is the position of
atom $i$ at time $t$. **Note** that fitting does not have to use the
same atoms as the calculation of the $RMSD$; *e.g.* a protein is
usually fitted on the backbone atoms (N, C$_{\alpha}$, C), but the
$RMSD$ can be computed of the backbone or of the whole protein.

Instead of comparing the structures to the initial structure at time
$t=0$ (so for example a crystal structure), one can also calculate
`eqn. %s <eqnrmsd>` with a structure at time $t_2=t_1-\tau$. This
gives some insight in the mobility as a function of $\tau$. A matrix
can also be made with the $RMSD$ as a function of $t_1$ and $t_2$,
which gives a nice graphical interpretation of a trajectory. If there
are transitions in a trajectory, they will clearly show up in such a
matrix.

Alternatively the $RMSD$ can be computed using a fit-free method with
the program `gmx rmsdist `:

$$
RMSD(t) ~=~     \left[\frac{1}{N^2}\sum_{i=1}^N \sum_{j=1}^N    \|{\bf r}_{ij}(t)-{\bf r}_{ij}(0)\|^2\right]^{\frac{1}{2}}
$$

where the *distance* **r**$_{ij}$ between atoms at time $t$ is
compared with the distance between the same atoms at time $0$.
