# gmxapi Python module reference

Version `gmxapi-version`.

The GROMACS Python package includes a
high-level scripting interface implemented in pure Python and a
lower-level API implemented as a C++ extension module. The pure Python
implementation provides the basic `gmxapi` module and classes with a
very stable syntax that can be maintained with maximal compatibility
while mapping to lower level interfaces that may take a while to sort
out. The separation also serves as a reminder that different execution
contexts may be implemented quite diffently, though Python scripts using
only the high-level interface should execute on all.

Package documentation is extracted from the `gmxapi` Python module and
is also available directly, using either `pydoc` from the command line
or `help` from within Python, such as during an interactive session.

Refer to the Python source code itself for additional clarification.

> **See also:** `gmxapi_package_documentation`


## Interface concepts

*gmxapi* commands return *references* to operations. Generally, the
operations are collected into a graph of data flow dependencies, and
only executed when the results are requested.

OperationReference


*gmxapi* uses a `Future` to reference an
operation output or data that may not yet be available.

Future


An `OperationReference` may provide several
named Futures on its *output* attribute.

A `Future` may be provided directly as
inputs to other *gmxapi* commands. *gmxapi* will execute the required
operation to get the data when it is needed.

To get an actual result in your Python script, you can call
`~Future.result()` on any *gmxapi* data reference. If the operation has
not yet been executed, the operation (and any operation dependencies)
will be executed immediately.

You can also force an operation to run by calling its
`~OperationReference.run()` method. But this is not generally necessary
unless your only goal is to produce output files on disk that are not
consumed in the same script.

In some cases, a `Future` can be
subscripted to get a new Future representing a slice of the original.
For instance, `commandline_operation`
objects have a *file* output that produces a mapping of command line
flags to output files (per the *output_files* parameter). This *file*
output can be subscripted with a single command line option to get a
`Future` for just one output file type. See
`gmxapi simulation preparation` for an illustrative example.

### Ensemble data flow

*gmxapi* automatically generates arrays of operations and parallel data
flow, when parallel inputs are provided to *gmxapi* command parameters.

When a `Future` represents the output of an
ensemble operation, `~Future.result()` returns a list with elements
corresponding to the ensemble members.

It is not currently possible to get a
`Future` for a specific ensemble member.

See `gmxapi ensemble` for more information.

## gmxapi basic package

    import gmxapi as gmx

gmxapi


function_wrapper


commandline_operation


subgraph


while_loop


## Simulation module

gmxapi.simulation


### Preparing simulations

read_tpr


gmxapi.simulation.read_tpr.OutputDataProxy


modify_input


gmxapi.simulation.modify_input.OutputDataProxy


### Running simulations

mdrun


gmxapi.simulation.mdrun.OutputDataProxy


## Utilities

gmxapi.utility


config


join_path


concatenate_lists


join_arrays


logical_not


make_constant


### Run time details

> [!NOTE]
> The *gmxapi.runtime* Python module is evolving. Some details are not
> yet well specified.

gmxapi.runtime


filtered_mpi_environ


filtered_prefixes


## Status messages and Logging

gmxapi.\_logging


## Exceptions module

gmxapi.exceptions


## gmx.version module

gmxapi.version


## Core API

gmxapi.\_gmxapi


### Exceptions

#### Module Exceptions

gmxapi.\_gmxapi.Exception

Root exception for the C++ extension module. Derives from
`gmxapi.exceptions.Error`.


FeatureNotAvailable


#### Wrapped C++ exceptions emitted through the supporting GROMACS library

gmxapi.\_gmxapi.MissingImplementationError


gmxapi.\_gmxapi.ProtocolError


gmxapi.\_gmxapi.UsageError


#### Other

No other C++ exceptions are expected, but will be wrapped in a
`Exception` to help tracing and reporting bugs.

gmxapi.\_gmxapi.UnknownException


### Functions

This documentation is provided for completeness and as an aid to
developers. Users of the `gmxapi` package, generally, should not need to
use the following tools directly.

#### Tools for launching simulations

gmxapi.\_gmxapi.from_tpr


gmxapi.\_gmxapi.create_context


#### Tools to manipulate TPR input files

gmxapi.\_gmxapi.copy_tprfile


gmxapi.\_gmxapi.read_tprfile


gmxapi.\_gmxapi.write_tprfile


gmxapi.\_gmxapi.rewrite_tprfile


#### Utilities

gmxapi.\_gmxapi.has_feature

Available features may depend on the package version, the details of the
supporting GROMACS installation, the software
environment detected when the package was built, or possibly on detected
runtime details. These feature checks are largely for internal use. The
`gmxapi` commands may adjust their behavior slightly depending on
feature checks, and (at worst) should produce meaningful error messages
or exceptions.

Named features:

- *create_context*: `create_context` can be
  used to initialize a `Context` with
  assigned resources.
- *mpi_bindings*: C++ extension module was built with `mpi4py`
  compatibility.


### Classes

gmxapi.\_gmxapi.Context


gmxapi.\_gmxapi.MDArgs


gmxapi.\_gmxapi.MDSession


gmxapi.\_gmxapi.MDSystem


gmxapi.\_gmxapi.SimulationParameters


gmxapi.\_gmxapi.TprFile


