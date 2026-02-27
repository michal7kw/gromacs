# gmxapi Python module reference

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

> **See also:** This copy of the documentation was built without
`installing the gmxapi package <../userguide/install>`, and therefore
lacks the full module reference. Refer to `gmxapi_package_documentation`
for instructions on building complete documentation, or [view
online](http://manual.gromacs.org/current/gmxapi/).

