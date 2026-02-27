# Public C++ API

## Overview

Trajectory analysis tools and pluggable MD extensions (such as code
based on the
[sample_restraint](https://gitlab.com/gromacs/gromacs/-/tree/main/python_packaging/sample_restraint)
example) use `gromacs/` headers supported by `libgromacs`.

html

For more on using `libgromacs` and the installed `gromacs/` headers,
refer to the [Legacy API](../doxygen/html-user/page_usinglibrary.xhtml)
documentation.

See also the [Trajectory Analysis
Framework](../doxygen/html-user/page_analysisframework.xhtml).


Software that uses new public API facilities (such as `gmxapi`) uses
CMake and `find_package(gmxapi)` to configure a build system to use the
`gmxapi/` headers and link to the library supporting the `::gmxapi` C++
namespace.

Currently, the `gmxapi` library conveys an indirect dependency on
`libgromacs`. Due to [a bug in CMake
3.24.0](https://gitlab.kitware.com/cmake/cmake/-/issues/23838),
`find_package(gmxapi)` must implicitly call
`find_package(gromacs${GROMACS_SUFFIX})` to avoid a spurious error, even
though client software does not generally need to explicitly use
`Gromacs::libgromacs` or its details.

## Client build system support

GROMACS relies heavily on
[CMake](https://cmake.org/documentation/) to configure and manage the
build system. The GROMACS installation
directly supports CMake configured client software through configuration
and "hints" files installed to `$GROMACS_ROOT/share/cmake/`.

`gmx --version` (or the appropriate `gmx$GROMACS_SUFFIX`) includes notes
on the original build toolchain that may or may not be sufficient for
configuring the client software build system.

### Compiler toolchain

Though not explicitly required, it is highly recommended that client
software build with a toolchain that closely matches that of the
GROMACS build to avoid binary
incompatibilities.

Each GROMACS installation (since 2022)
provides a CMake "hints" file that can be used to initialize your
`cmake` cache with the `-C`
[option](https://cmake.org/cmake/help/v3.24/manual/cmake.1.html#options).

For a GROMACS installation in
`$GROMACS_ROOT/`, the hints file for a given `GROMACS_SUFFIX` can be
found at
`$GROMACS_ROOT/share/cmake/gromacs$GROMACS_SUFFIX/gromacs-hints$GROMACS_SUFFIX.cmake`

The hints file is completely separate from the CMake configuration files
that support `find_package(gmxapi)` and `find_package(GROMACS)`.

However, using `-C path/to/gromacs-hints$GROMACS_SUFFIX.cmake` in your
client `cmake` configuration command line can help set appropriate
compiler options so that you have a better chance of building a
compatible binary. (I.e. it helps `gromacs_check_compiler` succeed.)

In addition to hints variables for CMake `enable_language` and
`find_package`, the hints file sets `GMX_CMAKE_VERSION` in case the
client build system needs to know the version of CMake that was used to
build the GROMACS installation.

### MPI support

GROMACS uses
[FindMPI](https://cmake.org/cmake/help/latest/module/FindMPI.html) (the
module that supports CMake `find_package(MPI ...)`) to locate and
configure compiler and linker options for MPI support. Client software
is advised to do the same.

If software support for MPI was detected by
GROMACS when built, the *gromacs-hints* file
(see above) will define input variables to help `find_package` locate
the same MPI installation.

### Caveats

If GROMACS is installed from a package built
in a different environment, the embedded toolchain information may be
inaccurate. This could make the `gmx --version` output misleading and
the *gromacs-hints* file useless. You may encounter spurious warnings
when configuring the client build system, and the client software may or
may not interact properly with the GROMACS
installation.

In a computing environment with multiple toolchains available (such as a
typical High Performance Computing (HPC) cluster), the toolchain may
depend on environment variables for consistent behavior. If environment
modules were used when setting up the GROMACS
build environment (e.g. `module load gcc openmpi/gcc`), it may be
necessary to load the same environment modules before building the
client software.

## `gmxapi` CMake package

The CMake configuration files installed with
GROMACS support the "Config mode" of CMake
[find_package](https://cmake.org/cmake/help/latest/command/find_package.html).
Unlike the `gromacs$GROMACS_SUFFIX` packages, CMake configuration files
only support a single `gmxapi` package name.

The `gmxapi` API and ABI hide most of the differences possible in
`libgromacs` from different build options. However, the
`gmxapi/mpi/resourceassignment.h` interface is affected by the original
choice of `GMX_MPI`. A stable interface is available to MPI-enabled
client software through the `gmxapi/mpi/gmxapi_mpi.h` template header.

Some GROMACS installations include multiple
builds. For instance, there may be a `libgromacs.so`, `libgromacs_d.so`,
`libgromacs_mpi.so`, and `libgromacs_mpi_d.so`, (according to build-time
values of `GMX_DOUBLE` and `GMX_MPI`) any *one* of which might be
provided by the `Gromacs::libgromacs` CMake target. Until resolution of
`4334`, only one version of the `Gromacs::gmxapi` is importable from a
GROMACS installation. Each
GROMACS installation (with `GMXAPI` `ON`)
overwrites the CMake configuration files for the previously installed
gmxapi support.

### Imported target

```cmake
Gromacs::gmxapi

The `gmxapi` package provides a single `Gromacs::gmxapi` target that
conveys access to the installed `gmxapi/` headers. The associated shared
object library will be differently named, depending on the build system
configuration options. (See `GMX_DOUBLE` and `GMX_MPI`).
```


## `gromacs` (and `gromacs$GROMACS_SUFFIX` packages)

The CMake machinery to support `find_package(GROMACS)` has two parts: a
`FindGROMACS.cmake` find module (found in
`share/gromacs/template/cmake/` in the installation and
`share/template/cmake/` in the source tree), and actual package
configuration files (`gromacs-config.cmake` and supporting files
installed to `share/cmake/` from input files in `src/gromacs/`).

`FindGROMACS.cmake` is a simple wrapper over the package configuration
files, providing a somewhat more convenient interface to the machinery
that supports multiple suffixed GROMACS
installations in the same installation prefix (see `GROMACS_SUFFIX`
variable below). This file is intended to be version-agnostic and remain
both forward-and backward-compatible even between major
GROMACS releases. All version-specific
information and the actual details about the compilation and linking
settings is in the package configuration files. Build systems willing to
utilize `FindGROMACS.cmake` can create a local copy of it and use it
like it is used in the installed
`share/gromacs/template/CMakeLists.txt`. The package configuration files
can also be used directly if desired, bypassing `FindGROMACS.cmake`.

When using `FindGROMACS.cmake`, `find_package(GROMACS)` is able to find
configurations for any of the `gromacs`, `gromacs_d`, `gromacs_mpi`, or
`gromacs_mpi_d` CMake package names. Otherwise, you must use the exact
package name that you are looking for. E.g. `find_package(gromacs_d)`.

### Imported targets

```cmake
Gromacs::libgromacs

Provides access to the installed core GROMACS
library and `gromacs/` headers:
`target_link_libraries(foo PRIVATE Gromacs::libgromacs)`.
```


```cmake
Gromacs::gmx

Represents the command line executable. For example, to set a local
CMake variable `_gmx_executable` to the executable path (with the
correct `GROMACS_SUFFIX`) you can use
`get_target_property(_gmx_executable Gromacs::gmx LOCATION)` in your
`CMakeLists.txt`
```


### Input options

Input options for influencing what to find

```cmake
GROMACS_SUFFIX

(only for `FindGROMACS.cmake`)

This CMake variable can be set before calling `find_package(GROMACS)` to
specify the GROMACS suffix to search for. If
not set, an unsuffixed version is searched for. If using the package
configuration files directly, the suffix must be set using
`find_package(GROMACS NAMES gromacs<suffix>)`.
```


```cmake
GROMACS_PREFER_STATIC

This CMake variable can be set before calling `find_package(GROMACS)` to
specify whether static or shared libraries are preferred if both are
available. It does not affect which GROMACS
installation is chosen, but if that installation has both static and
shared libraries available (installed from two different builds with the
same suffix), then this chooses the library to be returned in
`GROMACS_LIBRARIES`.
```


```cmake
GROMACS_DIR

This CMake (cache) variable is a standard mechanism provided by
`find_package`, and can be used to specify a hint where to search for
GROMACS. Also `CMAKE_PREFIX_PATH` can be used
for this purpose; see CMake documentation for `find_package` for more
details. `GROMACS_DIR` can also be set as an environment variable, and
this is done by `GMXRC`.
```


### Output variables

Output variables that specify how the found `libgromacs` and header
should be used:

```cmake
GROMACS_INCLUDE_DIRS

List of include directories necessary to compile against the
GROMACS headers. Currently, this includes the
path to GROMACS headers.
```


```cmake
GROMACS_LIBRARIES

List of libraries to link with to link against
GROMACS. Under the hood, this uses imported
CMake targets to represent `libgromacs`.
```


```cmake
GROMACS_DEFINITIONS

List of compile definitions (with `-D` in front) that are required to
compile the GROMACS headers.
```


```cmake
GROMACS_IS_DOUBLE

Whether the found GROMACS was compiled in
double precision.
```


```cmake
GROMACS_CXX_FLAGS

Required compiler flags.
```


### Macros/functions

Declared macros/functions that can be used for checking for correctness
of some settings:

```cmake
gromacs_check_double(GMX_DOUBLE)

Checks that the found GROMACS is in the
expected precision. The parameter `GMX_DOUBLE` should be the name of a
cache variable that specified whether double-precision was requested.
```


```cmake
gromacs_check_compiler(LANG)

Checks that the found GROMACS was compiled
with the same compiler that is used by the current CMake system.
Currently only `LANG=CXX` is supported.
```

