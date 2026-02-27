# C++ API

Public C++ application programming interfaces are available for
GROMACS installations depending on the
detected environment and user options when the
GROMACS build is configured with CMake.

- `gmxlibs`  
  - CMake target `Gromacs::gmxapi`, enabled by `GMXAPI` (default, when
    `BUILD_SHARED_LIBS` on non-Windows platforms), provides `gmxapi/`
    headers and `::gmxapi` C++ namespace.
  - CMake target `Gromacs::libgromacs`, enabled by
    `GMX_INSTALL_LEGACY_API` (default `OFF`), provides `gromacs/`
    headers and `::gmx` C++ namespace.

- `nblib`: Enabled by `GMX_INSTALL_NBLIB_API`. (default, when
  `BUILD_SHARED_LIBS` on non-Windows platforms)

html

- [Legacy API](../doxygen/html-user/index.xhtml): Enabled by
  `GMX_INSTALL_LEGACY_API`. (default `OFF`)


