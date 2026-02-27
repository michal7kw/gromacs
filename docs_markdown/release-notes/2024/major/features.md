# New and improved features

## The AWH exponential histogram growth can now be controlled

The accelerated weight histogram growth factor during the initial phase
was hard-coded to 3. Now this value can be controlled by the user. It is
set to 2 by default for increased stability.

If the TPR was generated with an earlier
GROMACS version, the old default value of 3
will be used.

## Added support for instrumentation based on wallcycle regions using NVTX/ROCTX/ITT

Basic support has been added for GPU tracing libraries so wallcycle main
and sub-regions will show up in tracing timelines which can help with
performance analysis. The tracing instrumentation support can be enabled
with one of the following CMake variables: `GMX_USE_NVTX`,
`GMX_USE_ROCTX`, `GMX_USE_ITT`.

`4446`

## Collective variables (Colvars) module support

The collective variables ([Colvars](https://colvars.github.io)) library
for enhanced sampling simulations has a new and improved interface,
which simplifies greatly its use and distribution with current and
future GROMACS releases. The new interface
requires *no patching* and supports a full integration of the Colvars
input and of its restart data with the GROMACS
TPR and CPT files, respectively.

For documentation and details, please see `this section ` of
the GROMACS doc along with the [Colvars doc
page](https://colvars.github.io/gromacs-2024/colvars-refman-gromacs.html)
for GROMACS. Additionally, messages in the
GROMACS discussion forum can also be tagged
with the [colvars keyword](https://gromacs.bioexcel.eu/tag/colvars) for
easier consultation.

## Automatic metric scaled AWH target distribution

The AWH target distribution can now be automatically scaled by sqrt(AWH
friction metric). Regions with higher friction (slower diffusion) will
get a higher target distribution. This should generally lower the
statistical error of the estimated free energy landscape. The new option
is called 'awh1-target-metric-scaling' and can be applied to further
modify all AWH target distributions and/or AWH user input, but is not
recommended in general in combination with Boltzmann or Local-Boltzmann
target distributions, due to the risk of feedback loops between the two
adaptive update mechanisms.
