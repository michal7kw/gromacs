# Curve fitting in GROMACS

## Sum of exponential functions

Sometimes it is useful to fit a curve to an analytical function, for
example in the case of autocorrelation functions with noisy tails.
GROMACS is not a general purpose curve-fitting
tool however and therefore GROMACS only
supports a limited number of functions. `Table %s `  lists
the available options with the corresponding command-line options. The
underlying routines for fitting use the Levenberg-Marquardt algorithm as
implemented in the lmfit package `162 <reflmfit>` (a bare-bones version
of which is included in GROMACS in which an
option for error-weighted fitting was implemented).

| Command line option | Functional form $f(t)$ | Note |
|----|----|----|
| exp | $e^{-t/{a_0}}$ |  |
| aexp | $a_1e^{-t/{a_0}}$ |  |
| exp_exp | $a_1e^{-t/{a_0}}+(1-a_1)e^{-t/{a_2}}$ | $a_2\ge a_0\ge 0$ |
| exp5 | $a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_4$ | $a_2\ge a_0\ge 0$ |
| exp7 | $a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_6$ | $a_4\ge a_2\ge a_0 \ge0$ |
| exp9 | $a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_7e^{-t/{a_6}}+a_8$ | $a_6\ge a_4\ge a_2\ge a_0\ge 0$ |

Overview of fitting functions supported in (most) analysis tools that
compute autocorrelation functions. The **Note** column describes
properties of the output parameters.

## Error estimation

Under the hood GROMACS implements some more
fitting functions, namely a function to estimate the error in
time-correlated data due to Hess `149 <refHess2002a>`:

$$
\varepsilon^2(t) =
2 \alpha\tau_1\left(1+\frac{\tau_1}{t}\left(e^{-t/\tau_1}-1\right)\right)
+ 2 (1-\alpha)\tau_2\left(1+\frac{\tau_2}{t}\left(e^{-t/\tau_2}-1\right)\right)
$$

where $\tau_1$ and $\tau_2$ are time constants (with
$\tau_2 \ge \tau_1$) and $\alpha$ usually is close to 1 (in the
fitting procedure it is enforced that $0\leq\alpha\leq 1$). This is
used in `gmx analyze ` for error estimation using

$$
\lim_{t\rightarrow\infty}\varepsilon(t) = \sigma\sqrt{\frac{2(\alpha\tau_1+(1-\alpha)\tau_2)}{T}}
$$

where $\sigma$ is the standard deviation of the data set and $T$ is
the total simulation time `149 <refHess2002a>`.

## Interphase boundary demarcation

In order to determine the position and width of an interface,
Steen-Sæthre *et al.* fitted a density profile to the following function

$$
f(x) ~=~ \frac{a_0+a_1}{2} - \frac{a_0-a_1}{2}{\rm
erf}\left(\frac{x-a_2}{a_3^2}\right)
$$

where $a_0$ and $a_1$ are densities of different phases, $x$ is
the coordinate normal to the interface, $a_2$ is the position of the
interface and $a_3$ is the width of the
interface `163 <refSteen-Saethre2014a>`. This is implemented in
`gmx densorder `.

## Transverse current autocorrelation function

In order to establish the transverse current autocorrelation function
(useful for computing viscosity  `164 <refPalmer1994a>`) the following
function is fitted:

$$
f(x) ~=~ e^{-\nu}\left({\rm cosh}(\omega\nu)+\frac{{\rm
sinh}(\omega\nu)}{\omega}\right)
$$

with $\nu = x/(2a_0)$ and $\omega = \sqrt{1-a_1}$. This is
implemented in `gmx tcaf `.

## Viscosity estimation from pressure autocorrelation function

The viscosity is a notoriously difficult property to extract from
simulations `149 <refHess2002a>`, `165 <refWensink2003a>`. It is *in
principle* possible to determine it by integrating the pressure
autocorrelation function `160 <refPSmith93c>`, however this is often
hampered by the noisy tail of the ACF. A workaround to this is fitting
the ACF to the following function `166 <refGuo2002b>`:

$$
f(t)/f(0) = (1-C) {\rm cos}(\omega t) e^{-(t/\tau_f)^{\beta_f}} + C
e^{-(t/\tau_s)^{\beta_s}}
$$

where $\omega$ is the frequency of rapid pressure oscillations (mainly
due to bonded forces in molecular simulations), $\tau_f$ and
$\beta_f$ are the time constant and exponent of fast relaxation in a
stretched-exponential approximation, $\tau_s$ and $\beta_s$ are
constants for slow relaxation and $C$ is the pre-factor that
determines the weight between fast and slow relaxation. After a fit, the
integral of the function $f(t)$ is used to compute the viscosity:

$$
\eta = \frac{V}{k_B T}\int_0^{\infty} f(t) dt
$$

This equation has been applied to computing the bulk and shear viscosity
using different elements from the pressure
tensor `167 <refFanourgakis2012a>`.
