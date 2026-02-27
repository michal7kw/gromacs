# Performance improvements

## Free-energy kernels on GPU with CUDA for perturbed non-bonded interactions

Free-energy calculations can now run faster on hardware with CUDA GPUs.
The benefits are largest when the CPU is relatively weak and/or there
are many perturbed non-bonded interactions. This feature has been
contributed by "MetaX Integrated Circuits Technology". Note that this
feature is currently marked as experimental.

## Optimized PME force transfers for multi-GPU runs with GPU-aware MPI

When using a GPU build with an (external) GPU-aware MPI library, the
transfer of PME force data between GPUs has been optimized to use
non-blocking MPI calls, improving performance by up to around 5%.

## More permissive AWH histogram equilibration check

The AWH histogram equilibration check, which is recommended to be used
with multiple walkers, used a too strict tolerance of 20% maximum
deviation. This tolerance can now be set by the user and the default
value is 30%. This value provide significantly faster convergence when
using multiple walkers. Since it does not deteriorate the convergence
with a single walkers, the histogram equilibration check is now turned
on by default.
