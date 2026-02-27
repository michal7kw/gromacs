# Portability

## MPI 3.0 is now required, 2.0 no longer supported

When building with MPI library, GROMACS 2026
requires it to be conformant with MPI 3.0. This should not be a problem
as long as you are using MPI library from the past decade (e.g., OpenMPI
1.8+, MPICH 3.0+).

`4933`

## Full support for HIP as GPU backend

The HIP backend can now be used for offloading all GPU kernels to AMD
devices. Ported kernels include direct GPU-to-GPU communication and
GPU-only run path.

`4947`

## MKL can now be found more easily from non-oneAPI distributions

The CMake support for finding and using MKL is now more flexible,
including now being able to find and link MKL when it was installed from
Debian's repositories.

`5508`
