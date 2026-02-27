# Codebase overview

The root directory of the GROMACS repository
only contains `CMakeLists.txt` (the root file for the CMake build
system), a few files supporting the build system, and a few standard
informative files (`README` etc.). The `INSTALL` is generated for source
packages from `docs/install-guide/index.rst`.

All other content is in the following top-level directories:

`admin/`  
Contains various scripts for developer use, as well as configuration
files and scripts for some of the tools used.

`api/`  
Contains code for the installable `/api/index`.

`cmake/`  
Contains code fragments and find modules for CMake. Some content here is
copied and/or adapted from newer versions of CMake than the minimum
currently supported. Default suppression file for valgrind is also
included here. See `build-system` for details of the build system.

`docs/`  
Contains the build system logic and source code for all documentation,
both user-facing and developer-facing. Some of the documentation is
generated from the source code under `src/`; see `dev-doc-layout`. This
directory also contains some developer scripts that use the Doxygen
documentation for their operation.

`scripts/`  
Contains the templates for `GMXRC` script, some other installed scripts,
as well as installation rules for all these scripts.

`share/`  
Contains data files that will be installed under `share/`. These include
a template for writing C++ analysis tools, and data files used by
GROMACS.

`src/`  
Contains all source code. See `dev-source-layout`.

`tests/`  
Contains build system logic for some high-level tests. Currently, only
the regression test build system logic, while other tests are under
`src/`.

## Source code organization

The sample code for the Trajectory Analysis Framework is in
`share/template/`.

Code for the gmxapi Python package and the sample MD extension module is
in `python_packaging/`.

`api/` holds code for the installable `/api/index`, including the legacy
`gromacs` headers and full sources for `libgmxapi` and `/nblib/index`.

The rest of the source code is under the `src/` directory.

The following figure shows a high-level view of components of what gets
built from the source code under `src/` and how the code is organized.
Arrows indicate the direction of dependencies. The build system is
described in detail in `build-system`. With default options, the green
and white components are built as part of the default target. The gray
parts are for testing, and are by default only built as part of the
`tests` target, but if `GMX_DEVELOPER_BUILD` is `ON`, then these are
included in the default build target. See `testutils` for details of the
testing side.

```dot
dev_high_level_components

concentrate = yes node \[ shape=box, style=filled, width=2 \]

subgraph {  
rank = same externals \[ label="externalsnsrc/external/", group=common,
style=rounded \] gtest \[ label="Google Test &
Mocknsrc/external/googletest/", group=test style="rounded,filled",
fillcolor="0 0 0.9" \]

} subgraph { rank = same libgromacs \[ label="libgromacsnsrc/gromacs/",
group=gmx, fillcolor="0.33 0.3 1" \] } testutils \[
label="testutilsnsrc/testutils/", group=test style="rounded,filled",
fillcolor="0 0 0.9" \] subgraph { rank = same gmx \[
label="gmxnsrc/programs/", group=gmx, fillcolor="0.33 0.3 1" \] tests \[
label="test binariesnsrc/.../tests/", group=test style="rounded,filled",
fillcolor="0 0 0.9" \] template \[ label="analysis
templatenshare/template/", group=common fillcolor="0.33 0.3 1" \]

> gmx -\> template \[ style=invis, constraint=no \]

}

libgromacs -\> externals gmx -\> libgromacs testutils -\> externals
testutils -\> gtest testutils -\> libgromacs tests -\> gtest tests -\>
libgromacs tests -\> testutils template -\> libgromacs
```


Only a few files related to the build system are included at the root
level. All actual code is in subdirectories:

`src/gromacs/`  
The code under this directory is built into a single library,
`libgromacs`. Installed headers are also located in this hierarchy. This
is the main part of the code, and is organized into further
subdirectories as *modules*. See below for details.

`src/programs/`  
The GROMACS executable `gmx` is built from
code under this directory. Also found here is some of the driver code
for the `mdrun` module called by `gmx`, and numerous end-to-end tests of
`gmx mdrun`.

`src/{...}/tests/`  
Various subdirectories under `src/` contain a subdirectory named
`tests/`. The code from each such directory is built into a test binary.
Some such directories also provide shared test code as object libraries
that is linked into multiple test binaries from different folders. See
`testutils` for details.

`src/testutils/`  
Contains shared utility code for writing Google Test tests. See
`testutils` for details.

`src/external/`  
Contains bundled source code for various libraries and components that
GROMACS uses internally. All the code from
these directories are built using our custom build rules into
`libgromacs`, or in some cases into the test binaries. Some CMake
options change which parts of this code are included in the build. See
`build-system` for some explanation about how the code in this directory
is used.

`src/external/build-fftw/`  
This folder contains the build system code for downloading and building
FFTW to be included into `libgromacs`.

When compiling, the include search path is set to `src/` by the
`legacy_modules` CMake target for many of the interfaces that do not
have clearer module ownership.

Some directories from under `src/external/` may also be included,
depending on the compilation options.

### Organization under `src/gromacs/`

The `libgromacs` library is built from code under `src/gromacs/`. Again,
the top-level directory contains build and installation rules for the
library.

The code is organized into subdirectories. These subdirectories are
denoted as `modules` throughout this documentation. Each module consists
of a set of routines that do some well-defined task or a collection of
tasks. Many modules are represented by distinct CMake targets, and
[target_link_libraries()](https://cmake.org/cmake/help/latest/command/target_link_libraries.html)
should be used to get access to the headers and linkable symbols. Other
modules are only expressed by the filesystem hierarchy, and their source
files are compiled directly into the monolithic `libgromacs` CMake
target.

Modules under `src/gromacs/` are not part of the public installed
interface. However, some of the headers that were traditionally
installed have been moved `api/legacy/include` (not duplicated in
`src/`) pending specification of an updated public API. These interfaces
are grouped into the `legacy_api` CMake target (in the build tree), and
are available through the IMPORTED `Gromacs::libgromacs` target for an
installation configured with `GMX_INSTALL_LEGACY_API=ON`. They are
installed into a corresponding hierarchy under `include/gromacs/` in the
installation directory.

> Historically, comments at the top of the header files contain a note
> about their visibility: public (installed), intra-library (can be used
> from inside the library), or intra-module/intra-file. Many of these
> comments remain, but they are no longer maintained.


All headers should compile by themselves, with installed headers doing
so without reference to variables defined in `config.h` or requiring
other headers to be included before it. No installed headers are allowed
to include `config.h`. Cyclic include dependencies prevent this, and
must be avoided because of this. This is best guaranteed by including
every header in some source file as the first header, even before
`config.h`.

Code inside the library should not unnecessarily include headers. In
particular, headers should not include other headers if a forward
declaration of a type is enough for the header. Within the library
source files, include only headers from other modules that are necessary
for that file. Check the `CMakeLists.txt` for the module to see whether
you need to
[target_link_libraries()](https://cmake.org/cmake/help/latest/command/target_link_libraries.html).
Many modules distinguish between a public interface and a private
interface intended only for use inside the module implementation. In
such cases, the public module headers (for use by other modules in the
library) are in `src/gromacs/modulename/include/gromacs/modulename`
subdirectories. Module private headers (located with the source files)
may be leaked into the include path, such as through the
`legacy_modules` target, but should not be used by other modules!

See `naming` for some common naming patterns for files that can help
locating declarations.

Tests, and data required for them, are in a `tests/` subdirectory under
the module directory. See `testutils` for more details.

## Documentation organization

All documentation (including this developer guide) is produced from
source files under `docs/`, except for some command-line help that is
generated from the source code (by executing the compiled `gmx` binary).
The build system provides various custom targets that build the
documentation; see `build-system` for details.

`docs/fragments/`  
Contains reStructuredText fragments used through `.. include::`
mechanism from various places in the documentation.

### User documentation

`docs/install-guide/`  
Contains reStructuredText source files for building the install guide
section of the user documentation, as well as the `INSTALL` file for the
source package. The build rules are in `docs/CMakeLists.txt`.

`docs/reference-manual/`  
Contains reStructuredText source files to generate the reference manual
for html and LaTeX.

`docs/manual/`  
Contains LaTeX helper files to build the reference (PDF) manual.

`docs/user-guide/`  
Contains reStructuredText source files used to build the user guide
section of the user documentation. The build rules are in
`docs/CMakeLists.txt`.

`docs/how-to/`  
Contains reStructuredText source files building the how-to section of
the user focused documentation.

### Unix man pages

Man pages for programs are generated by running the `gmx` executable
after compiling it, and then using Sphinx on the reStructuredText files
that `gmx` writes out.

The build rules for the man pages are in `docs/CMakeLists.txt`.

### Developer guide

`docs/dev-manual/`  
Contains reStructuredText source files used to build the developer
guide. The build rules are in `docs/CMakeLists.txt`.

The organization of the developer guide is explained on the
`front page of
the guide <dev guide>`.

### Doxygen documentation

`docs/doxygen/`  
Contains the build rules and some overview content for the Doxygen
documentation. See `doxygen` for details of how the Doxygen
documentation is built and organized.

> **TODO:** Create a separate page (at the front of the developer guide, and/or at
the main index.rst) that describes the documentation from readers'
perspective, and move relevant content there. This should contain just
an overview of how the documentation is organized in the source tree.


The Doxygen documentation is made of a few different parts. Use the list
below as a guideline on where to look for a particular kind of content.
Since the documentation has been written over a long period of time and
the approach has evolved, not all the documentation yet follows these
guidelines, but this is what we are aiming at.

documentation pages  
These contain mainly overview content, from general-level introduction
down into explanation of some particular areas of individual modules.
These are generally the place to start familiarizing with the code or a
new area of the code. They can be reached by links from the main page,
and also through cross-links from places in the documentation where that
information is relevant to understand the context.

module documentation  
These contain mainly technical content, explaining the general
implementation of a particular module and listing the classes, functions
etc. in the module. They complement pages that describe the concepts.
They can be reached from the Modules tab, and also from all individual
classes, functions etc. that make up the module.

class documentation  
These document the usage of an individual class, and in some cases that
of closely related classes. Where necessary (and time allowing), a
broader overview is given on a separate page and/or in the module
documentation.

method documentation  
These document the individual method. Typically, the class documentation
or other overview content is the place to look for how different methods
interact.

file and namespace documentation  
These are generally only placeholders for links, and do not contain much
else. The main content is the list of classes and other entities
declared in that file.
