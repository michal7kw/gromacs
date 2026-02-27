# Full installation instructions

Installation instructions for the `gmxapi` Python package, built on
GROMACS.

Command line examples assume the
[bash](https://www.gnu.org/software/bash/) shell.

> **Note:** Regarding multiple GROMACS installations

Many GROMACS users switch between multiple
GROMACS installations on the same computer
using an HPC module system and/or a
`GMXRC <getting access to |Gromacs|>` configuration script. For the
equivalent sort of environment switching with the `gmxapi` Python
package, we recommend installing it in a different [Python virtual
environment](https://www.google.com/search?q=python+virtual+environment)
for each GROMACS installation. Once built, a
particular copy of the `gmxapi` Python package always refers to the same
GROMACS installation.


> **Tip:** Unprivileged `pip install`

The following documentation contains frequent references to the
[pip](https://pip.pypa.io/en/stable/) tool for installing Python
packages. In some cases, an unprivileged user should use the `--user`
command line flag to tell [pip](https://pip.pypa.io/en/stable/) to
install packages into the user site-packages directory rather than the
default site-packages directory for the Python installation. This flag
is not appropriate when running `pip` in a virtual environment (as
recommended) and is omitted in this documentation. If you need the
`--user` flag, you should modify the example commands to look something
like `pip install --upgrade somepackage --user`


> **Info:** Python 3 executable names

These instructions use the executable names `python` and `pip` instead
of `python3` or `pip3`. Some Python installations require the `3`
suffix, but it is usually not necessary if you have already activated a
Python virtual environment (recommended).


## Overview

Typically, setting up the *gmxapi* Python package follows these three
steps. If this overview is sufficient for your computing environment,
you may disregard the rest of this document.

### Install GROMACS

Locate your GROMACS installation, or build and
install. GROMACS 2022 or higher is
recommended.

> **See also:** [GROMACS
installation](http://manual.gromacs.org/documentation/current/install-guide/index.html)


The following assumes GROMACS is installed to
`/path/to/gromacs`

### Set up a Python virtual environment

> **See also:** `gmxapi venv`


> [!NOTE]
> `mpi4py` may require additional arguments (compiler hints). See
> `mpi_requirements`

``` bash
python3 -m venv $HOME/myvenv
. $HOME/myvenv/bin/activate
python -m ensurepip --default-pip
pip install --upgrade pip setuptools wheel
pip install mpi4py
```

### Install the gmxapi Python package

Pull the `gmxapi` package from PyPI, build it for the
GROMACS installation at `/path/to/gromacs`,
and install it to the Python prefix for the current environment.

``` bash
. /path/to/gromacs/bin/GMXRC
pip install --no-cache-dir gmxapi
```

> **See also:** `installation`


## Background

*gmxapi* comes in three parts:

- GROMACS gmxapi library for C++.
- This Python package, supporting Python 3.9 and higher
- MD restraint plugins and sample gmxapi client code

### GROMACS requirements

The Python package requires a GROMACS
installation. Locate an existing GROMACS
installation, or [build and install
GROMACS](http://manual.gromacs.org/documentation/current/install-guide/index.html)
before proceeding.

> [!NOTE]
> Note that gmxapi requires that GROMACS is
> configured with `GMXAPI=ON` and `BUILD_SHARED_LIBS=ON`. These are
> enabled by default in most cases. If these options were overridden for
> your GROMACS installation, you will see
> CMake errors when trying to build and install the gmxapi Python
> package or other client software.

If your installation has a `GMXRC` file, "source" the file
`as you normally would <getting access to |Gromacs|>` before using
GROMACS. Otherwise, note the installation
location so that you can provide it when building the gmxapi package.

### Build system requirements

gmxapi can be built for Python 3.9 and higher.

You will need a C++ 17 compatible compiler and a reasonably up-to-date
version of CMake. Full gmxapi functionality may also require an MPI
compiler (e.g. `mpicc`).

Important: To build a module that can be imported by Python, you need a
Python installation that includes the Python headers. Unfortunately, it
is not always obvious whether these headers are present or where to find
them. The simplest answer is to just try to build the Python package
using these instructions, and if gmxapi is unable to find the Python
tools it needs, try a different Python installation or install the
additional development packages.

On a Linux system, this may require installing packages such as
`python-dev` and/or `python3-dev`. If you are building Python, either
from scratch or with a tool like `pyenv install` (see [wiki
entry](https://github.com/pyenv/pyenv/wiki#how-to-build-cpython-with---enable-shared)
), be sure to enable installation of the Python C library with the
`--enable-shared` flag. Alternatively, various Python distributions
provide a sufficient build environment while only requiring installation
into a user home directory. (Some examples below.)

If you are using an HPC system with software available through modules
you may be able to just `module load` a different Python installation
and find one that works.

### Python environment requirements

gmxapi requires Python 3.9 or higher. Check your version with
`python3 --version` or `python --version`.

> [!NOTE]
> The following documentation assumes you do not need to use a trailing
> '3' to access a Python 3 interpreter on your system. The default
> Python interpreter on your system may use `python3` and `pip3` instead
> of `python` and `pip`. You can check the version with
> `python3 --version` or `python --version` and `pip --version`.

To build and install, you need the Python packages for
[cmake](https://pypi.org/project/cmake/),
[networkx](https://pypi.org/project/networkx/), and
[setuptools](https://pypi.org/project/setuptools/) (all available from
[PyPI with pip](https://pip.pypa.io/en/stable/)).

For full functionality, you should also have
[mpi4py](https://pypi.org/project/mpi4py/) and
[numpy](https://www.numpy.org/). These requirements and version numbers
are listed in `requirements.txt`.

The easiest way to make sure you have the requirements installed, first
update [pip](https://pip.pypa.io/en/stable/), then use the
`requirements.txt` file provided with the repository. File paths in this
section are relative to the root directory of your local copy of the
GROMACS source.

Confirm that [pip](https://pip.pypa.io/en/stable/) is available, install
[pip](https://pip.pypa.io/en/stable/) if it is missing, or get
instructions on how to install [pip](https://pip.pypa.io/en/stable/):

``` bash
python -m ensurepip --default-pip
```

Install or upgrade required components:

``` bash
python -m pip install --upgrade pip
pip install --upgrade setuptools wheel
```

#### "requirements" files in GROMACS source tree

If you are building from source code in a local copy of the
GROMACS source repository, a
`requirements.txt` allows you to preinstall the Python requirements
before installing the `gmxapi` package.

> pip install -r python_packaging/gmxapi/requirements.txt

### Documentation build requirements

See `gmxapi_package_documentation`

### Testing requirements

Note that the test suite is only available in the
GROMACS source tree. (It is not part of the
installed package.) Acquire the GROMACS
sources with `git` or by downloading an archive, as documented
elsewhere.

Testing is performed with [pytest](https://docs.pytest.org/en/latest/).

`python_packaging/gmxapi/requirements.txt` lists additional requirements
for testing. With [pip](https://pip.pypa.io/en/stable/):

``` bash
pip install -r python_packaging/gmxapi/requirements.txt
```

To test the full functionality also requires an MPI parallel
environment. You will need the
[mpi4py](https://pypi.org/project/mpi4py/) Python package and an MPI
launcher (such as `mpiexec`, `mpirun`, a launcher provided by your HPC
queuing system, or whatever is provided by your favorite MPI package for
your operating system).

### MPI requirements

For the ensemble simulations features, you will need an MPI
installation.

On an HPC system, this means you will probably have to use `module load`
to load a compatible set of MPI tools and compilers. Check your HPC
documentation or try `module avail` to look for an `openmpi`, `mpich`,
or `mvapich` module and matching compiler module. This may be as simple
as:

``` bash
module load gcc
module load mpicc
```

If you are using a GROMACS installation that
is already available through `module load`, try to find a Python
installation with the `mpi4py` package that is also available through
`module load`. The *module* system will generally enforce toolchain
compatibility between the loaded modules. If you `module load` mpi4py or
a Python installation with mpi4py, you will probably want to use this
version of the package in your venv. (See `gmxapi venv`) If you
`module load` an MPI-enabled GROMACS
installation, `gmxapi` will try to check `mpi4py` for compatibility.

Note that the compilers loaded might not be the first compilers
discovered automatically by the build tools we will use below, so you
may have to specify compilers on the command line for consistency. It
may be necessary to require that GROMACS,
gmxapi, and the sample code are built with the same compiler(s).

Note that strange errors have been known to occur when
[mpi4py](https://pypi.org/project/mpi4py/) is built with a different
tool set than has been used to build Python and gmxapi. If the default
compilers on your system are not sufficient for
GROMACS or gmxapi, you may need to build,
e.g., OpenMPI or MPICH, and/or [build
mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html) with a
specific MPI compiler wrapper. This can complicate building in
environments such as [Conda](https://docs.conda.io/en/latest/). You
should be able to confirm that your MPI compiler wrapper is consistent
with your GROMACS tool chain by comparing the
output of `mpicc --version` with the compiler information reported by
`gmx --version`.

Set the `MPICC` environment variable to the MPI compiler wrapper and
forcibly reinstall [mpi4py](https://pypi.org/project/mpi4py/):

``` bash
export MPICC=`which mpicc`
pip install --no-cache-dir --upgrade --no-binary ":all:" --force-reinstall mpi4py
```

If you have a different MPI C compiler wrapper, substitute it for
`mpicc` above.

While `gmxapi` is configuring its build system during installation, it
will try to confirm the compatibility of the `mpi4py` toolchain with
that of the GROMACS installation. If they
appear incompatible, you should see a `CMake` message that includes a
guess at what you might try using for `MPICC`. (If using `pip`, consider
using the `--verbose` option for more build output.)

## Installing the Python package

We recommend using Python's [pip](https://pip.pypa.io/en/stable/)
package installer to automatically download, build, and install the
latest version of the gmxapi package into a Python [virtual
environment](https://docs.python.org/3/tutorial/venv.html), though it is
also possible to install without a virtual environment. If installing
without a virtual environment as an un-privileged user, you may need to
use the `--user` option with `pip install`.

### Recommended installation

The instructions in this section assume that *pip* is able to download
files from the internet. Alternatively, refer to
`gmxapi offline install`.

#### Locate or install GROMACS

You need a GROMACS installation that includes
the gmxapi headers and library.

> [!WARNING]
> gmxapi does not recognize multiple GROMACS
> installations to the same `CMAKE_INSTALL_PREFIX`.

> The Python package uses files installed to `.../share/cmake/gmxapi/`
> to configure its C++ component. These configuration files are
> overwritten when installing GROMACS to the
> same
> [CMAKE_INSTALL_PREFIX](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html).
> Overlapping GROMACS installations may occur
> when GROMACS is installed for multiple
> configurations of MPI support and floating point precision. (See
> `4334` and related issues.)

If GROMACS 2020 or higher is already
installed, *and* was configured with `GMXAPI=ON` at build time (the
default), you may be able to just source the
`GMXRC <getting access to |Gromacs|>` (so that the Python package knows
where to find GROMACS) and skip to the next
section. Note that some GROMACS installations,
such as in high-performance computing environments, may not install a
`GMXRC`, and may instead provide access to the
GROMACS installation through a
`module load gromacs` or similar command.

If necessary, install a supported version of
GROMACS. When building
GROMACS from source, be sure to configure
cmake with the flag `-DGMXAPI=ON` (default).

Set the environment variables for the GROMACS
installation so that the gmxapi headers and library can be found when
building the Python package. If you installed to a `gromacs-gmxapi`
directory in your home directory as above and you use the `bash` shell,
do:

``` bash
source $HOME/gromacs-gmxapi/bin/GMXRC
```

If you are using a GROMACS installation that
does not provide `GMXRC`, see [gmxapi cmake hints](#gmxapi cmake hints)
and additional CMake hints below.

#### Set up a Python virtual environment

We recommend installing the Python package in a virtual environment. If
not installing in a virtual environment, you may not be able to install
necessary prerequisites (e.g. if you are not an administrator of the
system you are on).

The following instructions use the `venv` module. Alternative virtual
environments, such as [Conda](https://docs.conda.io/en/latest/), should
work fine, but are beyond the scope of this document. (We welcome
contributed recipes!)

Depending on your computing environment, the Python 3 interpreter may be
accessed with the command `python` or `python3`. Use `python --version`
and `python3 --version` to figure out which you need to use. The
following assumes the Python 3 interpreter is accessed with `python3`.

> **Tip:** `--system-site-packages`  

It can be tricky to properly or optimally build MPI enabled software in
computing clusters, and administrators often provide prebuilt packages
like `mpi4py`. If your computing environment has multiple Python
installations, try to choose one that already includes `mpi4py`. When
you are using a Python installation that provides `mpi4py`, generally,
you should be sure to use the existing `mpi4py` installation in your new
virtual environment by creating the `venv` with the
`--system-site-packages` option.

In personal computing environments (laptops and workstations), it is
common to have multiple Python installations, and it can be hard to keep
packages in the different installations from conflicting with each
other. Unless you know that you want to inherit the `mpi4py` package
from the system installation, it is generally cleaner *not* to inherit
the system site-packages.

Create a Python 3 virtual environment:

``` bash
python3 -m venv $HOME/myvenv
```

*or* (see note):

``` bash
python3 -m venv --system-site-packages $HOME/myvenv
```

Activate the virtual environment. Your shell prompt will probably be
updated with the name of the environment you created to make it more
obvious.

``` none
$ source $HOME/myvenv/bin/activate
(myvenv)$
```

> [!NOTE]
> After activating the *venv*, `python` and `pip` are sufficient. (The
> '3' suffix will no longer be necessary and will be omitted in the rest
> of this document.)

Activating the virtual environment may change your shell prompt to
indicate the environment is active. The prompt is omitted from the
remaining examples, but the remaining examples assume the virtual
environment is still active. (Don't do it now, but you can deactivate
the environment by running `deactivate`.)

#### Install dependencies

It is always a good idea to update
[pip](https://pip.pypa.io/en/stable/),
[setuptools](https://pypi.org/project/setuptools/), and
[wheel](https://pypi.org/project/wheel/) before installing new Python
packages:

``` bash
pip install --upgrade pip setuptools wheel
```

The gmxapi installer requires a few additional packages. It is best to
make sure they are installed and up to date before proceeding.

``` bash
pip install --upgrade cmake pybind11
```

We use [mpi4py](https://pypi.org/project/mpi4py/) for some features and
to ensure compatible MPI bindings throughout your Python environment.
**If you did not inherit mpi4py from system site-packages** (see
`above <system-site-packages>`), make sure to [install
it](https://mpi4py.readthedocs.io/en/stable/install.html) using the same
MPI installation that we are building GROMACS
against, and build with compatible compilers.

``` bash
MPICC=`which mpicc` pip install --no-cache-dir --upgrade mpi4py
```

> **See also:** `mpi_requirements`


#### Install the latest version of gmxapi

Fetch and install the latest official version of gmxapi from the Python
Packaging Index. Avoid locally cached previously-built packages that may
be incompatible with your current environment or
GROMACS installation:

``` bash
# Get the latest official release.
pip install --no-cache-dir gmxapi
```

or:

``` bash
pip download gmxapi
pip install gmxapi-<version>.tar.gz
```

substituting the name of the downloaded source distribution archive.

> **Warning:** Avoid cached "wheel" packages.

`pip` downloads a source distribution archive for gmxapi, then builds a
"wheel" package for your GROMACS installation.
This "wheel" normally gets cached, and will be used by any later attempt
to `pip install gmxapi` instead of rebuilding. This is not what you
want, if you upgrade GROMACS or if you want to
install the Python package for a different
GROMACS configuration (e.g. double-precision
or different MPI option.)

You can use `--no-cache-dir` to force rebuild of the package and its
build dependencies. This may be slow, however, and you may want to use
cached dependencies. You can [avoid wheel
cache](https://pip.pypa.io/en/stable/topics/caching/#avoiding-caching)
for just one target package by installing from the filesystem instead of
directly from PyPI.

See also `4335`


The [PyPI repository](https://pypi.org/project/gmxapi/#history) may
include pre-release versions, but `pip` will ignore them unless you use
the `--pre` flag:

``` bash
# Get the latest version, including pre-release versions.
pip install --no-cache-dir --pre gmxapi
```

If `pip` does not find your GROMACS
installation, use one of the following environment variables to provide
a hint.

##### CMake hints

The `gmxapi` package is distributed with C++ source code that needs to
be compiled against GROMACS libraries. The
build system is configured using CMake, mediated by scikit-build-core.
Refer to [scikit-build-core
documentation](https://scikit-build-core.readthedocs.io/en/latest/configuration.html#configuring-cmake-arguments-and-defines)
for the best information on passing build options through the Python
package installer (*e.g.*\* [pip](https://pip.pypa.io/en/stable/)). (Be
sure to look at the "config-settings" and "Environment" tabs. The
"pyproject.toml" tabs is for package maintainers.)

**gmxapi_ROOT**

If you have a single GROMACS installation at
`/path/to/gromacs`, it is usually sufficient to provide this location to
`pip` through the `gmxapi_ROOT` environment variable, or as a CMake
variable definition.

environment

``` shell
gmxapi_ROOT=/path/to/gromacs pip install --no-cache-dir gmxapi
```


pip argument

``` shell
pip install --no-cache-dir gmxapi --config-setting=cmake.define.gmxapi_ROOT=/path/to/gromacs
```


**GROMACS CMake hints**

It can be important to use the same compiler tool chain for both
GROMACS and for client software (the Python
package C++ extension).

You can check `gmx --version` to see what compilers your installation
used, and make sure that you don't have incompatible compilers declared
by environment variables such as `$CXX`.

You can also use the GROMACS provided CMake
cache file to provide extra hints to the Python extension build system
about the software tools that were used to build
GROMACS. (For more information, read about the
`-C` [command line
option](https://cmake.org/cmake/help/latest/manual/cmake.1.html#options)
for CMake.)

In the following example,

- `${UNIQUE_PREFIX}` is the path to the directory that holds the
  GROMACS `bin`, `lib`, `share` directories,
  *etc*. It is *unique* because GROMACS
  provides CMake support for only one build configuration at a time
  through `.../share/cmake/gmxapi/`, even if there are multiple library
  configurations installed to the same location. See `4334`.
- `${SUFFIX}` is the suffix (e.g. *\_d*, *\_mpi*, etcetera) that
  distinguishes the particular build of
  GROMACS you want to target (refer to
  GROMACS installation instructions for more
  information.) `${SUFFIX}` may simply be empty, or `''`.

``` shell
pip install gmxapi \
  --config-settings=cmake.define.gmxapi_ROOT=${UNIQUE_PREFIX} \
  --config-settings=cmake.args=-C${UNIQUE_PREFIX}/share/cmake/gromacs${SUFFIX}/gromacs-hints${SUFFIX}.cmake
```

In sufficiently new `pip` versions, `-C` [is a shorter
alternative](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-C)
to `--config-settings=`. Do not confuse the `-C` option to `pip` with
the `-C` option to `cmake`.

> **See also:** [scikit-build-core
config-settings](https://scikit-build-core.readthedocs.io/en/latest/configuration.html#configuring-cmake-arguments-and-defines)


### Install from source

You can also install the `gmxapi` Python package from within a local
copy of the GROMACS source repository.
Assuming you have already obtained the GROMACS
source code and you are in the root directory of the source tree, you
will find the `gmxapi` Python package sources in the
`python_packaging/gmxapi` directory.

``` bash
cd python_packaging/gmxapi
pip install -r requirements.txt
pip install .
```

### Offline install

> **Tip:** Recommended, first:

`pip install --upgrade build pip setuptools wheel`


You can use `python -m build --skip-dependency-check` to build a binary
distribution archive (from the source distribution) for just the
*gmxapi* package, but then you will have to manually satisfy (separate)
dependencies in both the build and installation environments.

While you have internet access, you need to get access to the *gmxapi*
source distribution and the package dependencies. You will also want the
`wheel` and `build` packages in environments where the package(s) will
be built. Only `pip` is necessary once a gmxapi `wheel` is built.

The following instructions are paraphrased from
<https://pip.pypa.io/en/stable/user_guide/#installing-from-local-packages>

To build with internet access and then install without:

``` bash
# Remove any locally cached (previously built) wheels.
pip cache remove gmxapi

# Download gmxapi and dependencies from pypi.
pip wheel --wheel-dir DIR gmxapi
# or, using package source from the GROMACS repository
cd python_packaging/gmxapi
pip wheel --wheel-dir DIR .

# Later, install.
pip install --no-index --find-links=DIR DIR/gmxapi*whl
```

To download packages and dependencies for later build and installation:

``` bash
# if in the GROMACS source repository
cd python_packaging/gmxapi
# or download and expand the archive
pip download --destination-directory DIR gmxapi
tar xf DIR/gmxapi*
cd gmxapi*

# Pre-fetch dependencies to DIR
pip download --destination-directory DIR .

# Build and install from the source directory.
pip install --no-index --find-links=DIR .
```

### Building a source archive

A source archive for the gmxapi python package can be built from the
GROMACS source repository using the Python
[build](https://pypa-build.readthedocs.io/en/latest/) module.

Example:

``` bash
pip install --upgrade setuptools build
cd python_packaging/gmxapi
python -m build --sdist
```

This command will create a `dist` directory containing a source
distribution archive file. The file name has the form
`gmxapi-{version}.{suffix}`, where *version* is the version from the
package metadata, and *suffix* is an archive file extension determined
by the local environment and the current packaging specifications.

The version information is derived from `gmxapi.__version__` defined by
the `gmxapi.version` module. Pending refinement under `3851`, the gmxapi
version information is hard coded in the `version.py`. Make sure you
have an up-to-date version of the sources and that the version
information is appropriate before distributing a new release.

> **See also:** Python documentation for [creating a source
distribution](https://docs.python.org/3/distutils/sourcedist.html#creating-a-source-distribution)


Package maintainers may update the [online
repository](https://pypi.org/project/gmxapi/) by uploading a freshly
built `sdist` with
`python -m twine upload dist/gmxapi-{version}.{suffix}`. To update the
repository at the PyPI test server, use
`python -m twine upload --repository testpypi dist/gmxapi-{version}.{suffix}`.

## Accessing gmxapi documentation

Documentation for the Python classes and functions in the gmx module can
be accessed in the usual ways, using `pydoc` from the command line or
`help()` in an interactive Python session.

The complete documentation (which you are currently reading) can be
browsed [online](http://manual.gromacs.org/current/gmxapi/) or built
from a copy of the GROMACS source repository.

Documentation is built from a combination of Python module documentation
and static content, and requires a local copy of the
GROMACS source repository.

### Build with GROMACS

To build the full gmxapi documentation with
GROMACS, configure
GROMACS with `-DGMX_PYTHON_PACKAGE=ON` and
build the GROMACS documentation normally. This
will first build the *gmxapi* Python package and install it to a
temporary location in the build tree. Sphinx can then import the package
to automatically extract Python docstrings.

Note that this is an entirely CMake-driven installation and Python
dependencies will not be installed automatically. You can update your
Python environment (before configuring with CMake) using the
`requirements.txt` files provided in the `python_packaging/` directory
of the repository. Example:

``` bash
pip install -r python_packaging/gmxapi/requirements.txt
```

Sometimes the build environment can choose a different Python
interpreter than the one you intended. You can set the
`Python3_ROOT_DIR` or `CMAKE_PREFIX_PATH` CMake variable to explicitly
choose the Python installation or *venv* directory. See also [CMake
FindPython3](https://cmake.org/cmake/help/latest/module/FindPython3.html).

If you use pyenv or pyenv-virtualenv to dynamically manage your Python
version, you can help identify a particular version with
`pyenv version-name` and the directory with `pyenv prefix {version}`.
For example:

``` bash
-DPython3_ROOT_DIR=$(pyenv prefix $(pyenv version-name))
```

> **TODO:** Document sample_restraint package. Reference `3027`


## Testing

Note [testing requirements](#testing-requirements) above.

After installing the `gmxapi` Python package, you can run the Python
test suite from the GROMACS source tree.
Example:

``` bash
# Assuming you are in the root directory of the repository:
pytest python_packaging/gmxapi/test/
```

Refer to `python_packaging/README.md` for more detailed information.

## Troubleshooting

### ImportError at run time with dynamic linking error

Symptom: Python fails with a weird `ImportError` citing something like
`dlopen`:

``` bash
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ImportError: dlopen(/.../gmxapi/_gmxapi.so, 0x0002): Symbol not found:
__ZN12gmxapicompat11readTprFileERKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
  Referenced from: /.../gmxapi/_gmxapi.so
  Expected in: /path/to/gromacs/lib/libgmxapi_mpi_d.0.3.1.dylib
```

Inconsistencies in the build and run time environments can cause dynamic
linking problems at run time. This could occur if you reinstall
GROMACS built with a different compiler, or if
`pip` or `CMake` somehow get tricked into using the wrong compiler tool
chain.

Refer to the [gmxapi cmake hints](#gmxapi cmake hints) for notes about
compiler toolchains. Rebuild and reinstall the gmxapi Python package
with `--no-cache-dir` and provide the `gromacs-hints.cmake` file for the
GROMACS installation you intend to use.

### AttributeError: module 'enum' has no attribute 'IntFlag'

If you had older versions of some of the dependencies installed, you
might have picked up a transitive dependency on the `enum34` package.
Try:

``` bash
pip uninstall -y enum34
```

and see if that fixes the problem. If not, try a fresh virtual
environment (see above) to help narrow down the problem before you [open
an issue](https://gitlab.com/gromacs/gromacs/-/issues/).

### Errors regarding pybind11

An error may occur in `setup.py` with output that contains something
like the following:

``` bash
ModuleNotFoundError: No module named 'pybind11'
Building wheel for gmxapi (pyproject.toml): finished with status 'error'
ERROR: Failed building wheel for gmxapi
Failed to build gmxapi
ERROR: Could not build wheels for gmxapi, which is required to install pyproject.toml-based projects
```

The important information here is that `pybind11` was not found.

Build dependencies aren't always automatically installed. Even if you
are using `pip`, you may have disabled automatic dependency fulfillment
with an option like `--no-build-isolation` or `--no-deps`.

In any case, the problem should be resolved by explicitly installing the
`pybind11` Python package before attempting to build `gmxapi`:

``` bash
pip install --upgrade pybind11
```

### Couldn't find the `gmxapi` support library?

If you don't want to "source" your `GMXRC <getting access to |Gromacs|>`
file, you can tell the package where to find a gmxapi compatible
GROMACS installation with `gmxapi_ROOT`. E.g.
`gmxapi_ROOT=/path/to/gromacs pip install .`

Before updating the `gmxapi` package it is generally a good idea to
remove the previous installation and to start with a fresh build
directory. You should be able to just `pip uninstall gmxapi`.

Do you see something like the following?

``` none
CMake Error at gmx/core/CMakeLists.txt:45 (find_package):
   Could not find a package configuration file provided by "gmxapi" with any
   of the following names:

     gmxapiConfig.cmake
     gmxapi-config.cmake

   Add the installation prefix of "gmxapi" to CMAKE_PREFIX_PATH or set
   "gmxapi_ROOT" to a directory containing one of the above files.  If "gmxapi"
   provides a separate development package or SDK, be sure it has been
   installed.
```

This could be because

- GROMACS is not already installed
- GROMACS was built without the CMake variable
  `GMXAPI=ON`
- or if `gmxapi_ROOT` (or `GROMACS_DIR`) is not a path containing
  directories like `bin` and `share`.

If you are not a system administrator you are encouraged to install in a
Python virtual environment, created with virtualenv or
[Conda](https://docs.conda.io/en/latest/). Otherwise, you will need to
specify the `--user` flag to `pip`.

Two of the easiest problems to run into are incompatible compilers and
incompatible Python. Try to make sure that you use the same C and C++
compilers for GROMACS, for the Python package,
and for the sample plugin. These compilers should also correspond to the
`mpicc` compiler wrapper used to [compile
mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html). In order
to build the Python package, you will need the Python headers or
development installation, which might not already be installed on the
machine you are using. (If not, then you will get an error about missing
`Python.h` at some point.) If you have multiple Python installations (or
modules available on an HPC system), you could try one of the other
Python installations, or you or a system administrator could install an
appropriate Python dev package. Alternatively, you might try installing
your own Anaconda or MiniConda in your home directory.

If an attempted installation fails with CMake errors about missing
“gmxapi”, make sure that GROMACS is installed
and can be found during installation. For instance,

``` bash
gmxapi_ROOT=/Users/eric/gromacs pip install --verbose gmxapi
```

Pip and related Python package management tools can be a little too
flexible and ambiguous sometimes. If things get really messed up, try
explicitly uninstalling the `gmxapi` module and its dependencies, then
do it again and repeat until `pip` can no longer find any version of any
of the packages.

``` bash
pip uninstall gmxapi
pip uninstall cmake
# ...
```

Successfully running the test suite is not essential to having a working
`gmxapi` package. We are working to make the testing more robust, but
right now the test suite is a bit delicate and may not work right, even
though you have a successfully built the `gmxapi` package. If you want
to troubleshoot, though, the main problems seem to be that automatic
installation of required python packages may not work (requiring manual
installations, such as with `pip install somepackage`) and ambiguities
between python versions.

If you are working in a development branch of the repository, note that
the upstream branch may be reset to `main` after a new release is
tagged. In general, but particularly on the `devel` branch, when you do
a `git pull`, you should use the `--rebase` flag.

If you fetch this repository and then see a git status like this:

``` bash
$ git status
On branch devel
Your branch and 'origin/devel' have diverged,
and have 31 and 29 different commits each, respectively.
```

then `gmxapi` has probably entered a new development cycle. You can do
`git pull --rebase` to update to the latest development branch.

If you do a `git pull` while in `devel` and get a bunch of unexpected
merge conflicts, do `git merge --abort; git pull --rebase` and you
should be back on track.

If you are developing code for gmxapi, this should be an indication to
rebase your feature branches for the new development cycle.
