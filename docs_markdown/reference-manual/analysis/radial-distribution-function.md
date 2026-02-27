# Radial distribution functions

`gmx rdf `  
The *radial distribution function* (RDF) or pair correlation function
$g_{AB}(r)$ between particles of type $A$ and $B$ is defined in
the following way:

$$
\begin{aligned}
\begin{array}{rcl}
g_{AB}(r)&=&    {\displaystyle \frac{\langle \rho_B(r) \rangle}{\langle\rho_B\rangle_{local}}}         \\
&=&    {\displaystyle \frac{1}{\langle\rho_B\rangle_{local}}}{\displaystyle \frac{1}{N_A}}
\sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B}
{\displaystyle \frac{\delta( r_{ij} - r )}{4 \pi r^2}}         \\
\end{array}
\end{aligned}
$$

with $\langle\rho_B(r)\rangle$ the particle density of type $B$ at a
distance $r$ around particles $A$, and
$\langle\rho_B\rangle_{local}$ the particle density of type $B$
averaged over all spheres around particles $A$ with radius $r_{max}$
(see `Fig. %s <fig-rdfex>` C).

![plots/rdf.*](plots/rdf.*)

*Definition of slices in gmx rdf &lt;gmx rdf&gt;: A. $gAB(r)$. B. $gAB(r, θ)$. The slices are colored gray. C. Normalization $⟨ρB⟩local$. D. Normalization $⟨ρB⟩local, θ$. Normalization volumes are colored gray.*

Usually the value of $r_{max}$ is half of the box length. The
averaging is also performed in time. In practice the analysis program
`gmx rdf ` divides the system into spherical slices (from $r$
to $r+dr$, see `Fig. %s <fig-rdfex>` A) and makes a histogram in stead
of the $\delta$-function. An example of the RDF of oxygen-oxygen in
SPC water `80 <refBerendsen81>` is given in `Fig. %s <fig-rdf>`

![plots/rdfO-O.*](plots/rdfO-O.*)

*$gOO(r)$ for Oxygen-Oxygen of SPC-water.*

With `gmx rdf ` it is also possible to calculate an angle
dependent rdf $g_{AB}(r,\theta)$, where the angle $\theta$ is
defined with respect to a certain laboratory axis ${\bf e}$, see
`Fig. %s <fig-rdfex>` B.

$$
g_{AB}(r,\theta) = {1 \over \langle\rho_B\rangle_{local,\:\theta }}
{1 \over N_A} \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} {\delta( r_{ij} - r )
\delta(\theta_{ij} -\theta) \over 2 \pi r^2 sin(\theta)}
$$

$$
cos(\theta_{ij}) = {{\bf r}_{ij} \cdot {\bf e} \over \|r_{ij}\| \;\| e\| }
$$

This $g_{AB}(r,\theta)$ is useful for analyzing anisotropic systems.
**Note** that in this case the normalization
$\langle\rho_B\rangle_{local,\:\theta}$ is the average density in all
angle slices from $\theta$ to $\theta + d\theta$ up to $r_{max}$,
so angle dependent, see `Fig. %s <fig-rdfex>` D.
