# Portability

## Added support for Hygon Dhyana CPU architecture

Support for hardware detection and related heuristics has been
implemented for the Hygon Dhyana derived from the first-gen AMD Zen
which it shares most of its architectural details with.

## Enabled PME offload support with OpenCL on NVIDIA and Intel GPUs

Thanks to portability improvements, the previously disabled PME OpenCL
offload is now enabled also on NVIDIA and Intel GPUs.

## Fixed building on Solaris with GCC

GROMACS now builds on Solaris with GCC (tested
on illumos distribution openindiana, "Hipster" rolling release, using
GCC 5, 6, 7, and 8).
