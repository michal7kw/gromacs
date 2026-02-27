# Portability

## Always use the Boost version bundled with GROMACS

Boost 1.83 is known to have compatibility issues when using Clang
compiler on FreeBSD and Linux. This changes to always use the bundled
Boost version, even if another version is present on the system.

`4893`
