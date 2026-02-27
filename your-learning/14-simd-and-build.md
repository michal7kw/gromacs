# Module 14: SIMD Abstraction and Build System

> **Prerequisites:** [Module 6 (nbnxm)](06-nbnxm-nonbonded.md), [Module 10 (C++ Patterns)](10-cpp-patterns.md), [Module 13 (GPU)](13-gpu-acceleration.md)
> **Key files:** `src/gromacs/simd/include/gromacs/simd/simd.h`, `src/gromacs/simd/include/gromacs/simd/simd_math.h`, `CMakeLists.txt`, `cmake/gmxDetectSimd.cmake`, `cmake/gmxManageMPI.cmake`

---

## 14.1 Why SIMD Matters for MD

SIMD (Single Instruction, Multiple Data) processes 4-16 floating-point values simultaneously. For non-bonded force kernels that evaluate millions of atom pairs per step, this can provide 2-8x speedups:

```
Scalar:  [r1²] → [1/r] → [F1]      [r2²] → [1/r] → [F2]      (sequential)
SIMD-4:  [r1², r2², r3², r4²] → [1/r, 1/r, 1/r, 1/r] → [F1, F2, F3, F4]  (parallel)
```

GROMACS supports 10+ SIMD instruction sets across x86, ARM, and IBM Power. Rather than writing separate code for each, it uses a **portable abstraction layer**.

---

## 14.2 The SIMD Type System

**File:** `src/gromacs/simd/include/gromacs/simd/simd.h`

### Core Types

| Type | Purpose | Width (AVX2-256) |
|------|---------|-----------------|
| `SimdFloat` | Single-precision vector | 8 floats |
| `SimdDouble` | Double-precision vector | 4 doubles |
| `SimdFInt32` | Integer vector (float-width) | 8 ints |
| `SimdDInt32` | Integer vector (double-width) | 4 ints |
| `SimdFBool` | Boolean vector (float-width) | 8 bools |
| `SimdDBool` | Boolean vector (double-width) | 4 bools |

### Backend Selection (Compile Time)

```cpp
// simd.h (lines 128-152) — preprocessor dispatch
#if GMX_SIMD_X86_SSE2
#    include "impl_x86_sse2/impl_x86_sse2.h"
#elif GMX_SIMD_X86_SSE4_1
#    include "impl_x86_sse4_1/impl_x86_sse4_1.h"
#elif GMX_SIMD_X86_AVX_128_FMA
#    include "impl_x86_avx_128_fma/impl_x86_avx_128_fma.h"
#elif GMX_SIMD_X86_AVX_256
#    include "impl_x86_avx_256/impl_x86_avx_256.h"
#elif GMX_SIMD_X86_AVX2_256
#    include "impl_x86_avx2_256/impl_x86_avx2_256.h"
#elif GMX_SIMD_X86_AVX2_128
#    include "impl_x86_avx2_128/impl_x86_avx2_128.h"
#elif GMX_SIMD_X86_AVX_512
#    include "impl_x86_avx_512/impl_x86_avx_512.h"
#elif GMX_SIMD_ARM_NEON_ASIMD
#    include "impl_arm_neon_asimd/impl_arm_neon_asimd.h"
#elif GMX_SIMD_ARM_SVE
#    include "impl_arm_sve/impl_arm_sve.h"
#elif GMX_SIMD_IBM_VSX
#    include "impl_ibm_vsx/impl_ibm_vsx.h"
#elif GMX_SIMD_REFERENCE
#    include "impl_reference/impl_reference.h"
#else
#    include "impl_none/impl_none.h"
#endif
```

### Type Traits for Generic Code

```cpp
// simd.h (lines 424-485)
struct SimdFloatTag {};
struct SimdDoubleTag {};

template<> struct SimdTraits<SimdFloat> {
    using type = float;
    static constexpr int width = GMX_SIMD_FLOAT_WIDTH;
    using tag = SimdFloatTag;
};
```

---

## 14.3 How Backends Implement SimdFloat

Each backend defines `SimdFloat` differently but exposes the same operations:

### Reference Implementation (Pure C++)

**File:** `src/gromacs/simd/include/gromacs/simd/impl_reference/impl_reference_simd_float.h`

```cpp
class SimdFloat {
public:
    SimdFloat() {}
    SimdFloat(float f) { simdInternal_.fill(f); }
    std::array<float, GMX_SIMD_FLOAT_WIDTH> simdInternal_;
};

static inline SimdFloat gmx_simdcall simdLoad(const float* m, SimdFloatTag = {}) {
    SimdFloat a;
    std::copy(m, m + a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

static inline SimdFloat gmx_simdcall operator+(SimdFloat a, SimdFloat b) {
    SimdFloat res;
    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
        res.simdInternal_[i] = a.simdInternal_[i] + b.simdInternal_[i];
    return res;
}

static inline SimdFloat gmx_simdcall fma(SimdFloat a, SimdFloat b, SimdFloat c) {
    return a * b + c;  // Two operations (no hardware FMA)
}
```

### SSE2 Implementation (x86, 128-bit)

**File:** `src/gromacs/simd/include/gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd_float.h`

```cpp
class SimdFloat {
public:
    SimdFloat() {}
    SimdFloat(float f) : simdInternal_(_mm_set1_ps(f)) {}
    SimdFloat(__m128 simd) : simdInternal_(simd) {}
    __m128 simdInternal_;  // 4 floats
};

static inline SimdFloat gmx_simdcall simdLoad(const float* m, SimdFloatTag = {}) {
    assert(std::size_t(m) % 16 == 0);  // Must be 16-byte aligned
    return { _mm_load_ps(m) };
}

static inline SimdFloat gmx_simdcall operator+(SimdFloat a, SimdFloat b) {
    return { _mm_add_ps(a.simdInternal_, b.simdInternal_) };
}

// SSE2 has no FMA — emulated with multiply + add
static inline SimdFloat gmx_simdcall fma(SimdFloat a, SimdFloat b, SimdFloat c) {
    return { _mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_) };
}

// SSE2 has hardware rsqrt (low precision)
static inline SimdFloat gmx_simdcall rsqrt(SimdFloat x) {
    return { _mm_rsqrt_ps(x.simdInternal_) };
}
```

### AVX2-256 Implementation (x86, 256-bit)

**File:** `src/gromacs/simd/include/gromacs/simd/impl_x86_avx2_256/impl_x86_avx2_256_simd_float.h`

```cpp
// AVX2 has TRUE hardware FMA (single instruction, better precision)
static inline SimdFloat gmx_simdcall fma(SimdFloat a, SimdFloat b, SimdFloat c) {
    return { _mm256_fmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fms(SimdFloat a, SimdFloat b, SimdFloat c) {
    return { _mm256_fmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fnma(SimdFloat a, SimdFloat b, SimdFloat c) {
    return { _mm256_fnmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}
```

### Key Difference: FMA

| Backend | FMA | Instruction | Rounding |
|---------|-----|-------------|----------|
| Reference | Emulated | `a*b + c` (2 ops) | Double rounding |
| SSE2 | Emulated | `_mm_add_ps(_mm_mul_ps())` | Double rounding |
| AVX-128-FMA | Hardware | `_mm_fmadd_ps` | Single rounding |
| AVX2-256 | Hardware | `_mm256_fmadd_ps` | Single rounding |
| AVX-512 | Hardware | `_mm512_fmadd_ps` | Single rounding |
| ARM NEON | Hardware | `vfmaq_f32` | Single rounding |

Hardware FMA gives both better performance (1 instruction vs 2) and better precision (single rounding instead of double rounding).

---

## 14.4 SIMD Math Functions

**File:** `src/gromacs/simd/include/gromacs/simd/simd_math.h`

MD requires functions like `1/sqrt(r²)`, `erfc()`, `exp()`, and `log()` on SIMD vectors. These are implemented using Newton-Raphson iterations and polynomial approximations.

### Inverse Square Root (the most critical function)

```cpp
// Newton-Raphson iteration for rsqrt
static inline SimdFloat gmx_simdcall rsqrtIter(SimdFloat lu, SimdFloat x) {
    SimdFloat tmp1 = x * lu;
    SimdFloat tmp2 = SimdFloat(-0.5F) * lu;
    tmp1 = fma(tmp1, lu, SimdFloat(-3.0F));
    return tmp1 * tmp2;
}

// Full-precision invsqrt
static inline SimdFloat gmx_simdcall invsqrt(SimdFloat x) {
    SimdFloat lu = rsqrt(x);  // Hardware low-precision estimate
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);   // First Newton-Raphson iteration
#endif
#if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);   // Second iteration (if needed)
#endif
    return lu;
}
```

**How it works:**
1. Hardware `rsqrt` gives ~12 bits of precision (SSE) or ~14 bits (AVX-512)
2. Each Newton-Raphson iteration roughly doubles the precision
3. For single-precision (24 bits), SSE needs 2 iterations; AVX-512 needs 1

### Logarithm (Base-2) via Minimax Polynomial

```cpp
static inline SimdFloat gmx_simdcall log2(SimdFloat x) {
    // Coefficients from minimax polynomial fit
    const SimdFloat CL9(0.342149508897807708152F);
    const SimdFloat CL7(0.411570606888219447939F);
    const SimdFloat CL5(0.577085979152320294183F);
    const SimdFloat CL3(0.961796550607099898222F);
    const SimdFloat CL1(2.885390081777926774009F);

    SimdFloat fExp, x2, p;
    SimdFInt32 iExp;

    // Separate mantissa and exponent using frexp
    x = frexp(x, &iExp);
    fExp = cvtI2R(iExp);

    // Range reduction to [-1, 1]
    x = (x - one) * inv(x + one);
    x2 = x * x;

    // Horner's method polynomial evaluation
    p = fma(CL9, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, fExp);

    return p;
}
```

### Safe vs Unsafe Math

```cpp
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall sqrt(SimdFloat x) {
    if constexpr (opt == MathOptimization::Safe) {
        SimdFloat res = maskzInvsqrt(x, setZero() < x);  // Handles x=0
        return res * x;
    } else {
        return x * invsqrt(x);  // Faster, but x must be > 0
    }
}
```

The `Safe` mode adds a mask to handle zeros; `Unsafe` assumes valid input. Non-bonded kernels use `Unsafe` since `r²` is always positive for interacting pairs.

---

## 14.5 SIMD in Non-Bonded Kernels

The nbnxm kernel layouts (Module 6) are designed specifically for SIMD:

### 4xM Layout

```
i-cluster: 4 atoms (always 4)
j-cluster: M atoms (M = SIMD float width)

For AVX2-256 (M=8):
  Each atom in i interacts with 8 atoms in j
  → 4 × 8 = 32 interactions per cluster pair
  → 8 LJ+Coulomb evaluations happen in a single SIMD instruction
```

### 2xMM Layout

```
i-cluster: 2 atoms
j-cluster: 2M atoms (double-width)

For AVX2-256 (M=8, MM=16):
  Each atom in i interacts with 16 atoms in j
  → 2 × 16 = 32 interactions per cluster pair
```

### SIMD Kernel Pseudocode

```cpp
// Simplified inner loop of a 4xM non-bonded kernel
for (int jIndex = 0; jIndex < numJClusters; jIndex++) {
    // Load 4 atom positions from i-cluster (scalar, then broadcast)
    SimdFloat ix0 = SimdFloat(xi[0][XX]);
    SimdFloat iy0 = SimdFloat(xi[0][YY]);
    SimdFloat iz0 = SimdFloat(xi[0][ZZ]);

    // Load M atom positions from j-cluster (SIMD load)
    SimdFloat jx = simdLoad(xj + XX*stride);
    SimdFloat jy = simdLoad(xj + YY*stride);
    SimdFloat jz = simdLoad(xj + ZZ*stride);

    // Distance calculation (M distances at once)
    SimdFloat dx = ix0 - jx;
    SimdFloat dy = iy0 - jy;
    SimdFloat dz = iz0 - jz;
    SimdFloat rsq = fma(dx, dx, fma(dy, dy, dz * dz));

    // 1/r calculation (M values at once)
    SimdFloat rinv = invsqrt(rsq);
    SimdFloat rinvsq = rinv * rinv;

    // LJ force: 6*C6/r^8 - 12*C12/r^14
    SimdFloat rinvsix = rinvsq * rinvsq * rinvsq;
    SimdFloat FrLJ = fma(c12 * rinvsix, rinvsix, fneg(c6) * rinvsix);

    // Coulomb force
    SimdFloat FrCoul = qq * rinv * rinvsq;

    // Accumulate forces
    SimdFloat Ftotal = (FrLJ + FrCoul) * rinvsq;
    fix0 = fma(Ftotal, dx, fix0);  // M force updates at once
    fiy0 = fma(Ftotal, dy, fiy0);
    fiz0 = fma(Ftotal, dz, fiz0);
}
```

---

## 14.6 CMake Build System Overview

**File:** `CMakeLists.txt` (top-level)

### Major Build Options

```cmake
# SIMD instruction set (lines 241-280)
gmx_option_multichoice(GMX_SIMD
    "SIMD instruction set for CPU kernels"
    "AUTO"
    AUTO None SSE2 SSE4.1 AVX_128_FMA AVX_256 AVX2_256
    AVX2_128 AVX_512 ARM_NEON_ASIMD ARM_SVE IBM_VSX Reference)

# GPU framework (lines 242-274)
gmx_option_multichoice(GMX_GPU
    "Framework for GPU acceleration"
    OFF
    OFF CUDA OpenCL SYCL HIP)

# FFT library (lines 288-304)
gmx_option_multichoice(GMX_FFT_LIBRARY
    "FFT library"
    "${GMX_FFT_LIBRARY_DEFAULT}"
    fftw3 mkl "fftpack[built-in]")

# MPI (lines 220-221)
option(GMX_MPI "Build with real MPI" OFF)
option(GMX_THREAD_MPI "Build with thread-MPI" ON)

# C++ standard (lines 69-76)
set(CMAKE_CXX_STANDARD 17)
```

### Compiler Requirements

```cmake
# Minimum compiler versions (lines 86-92)
# Clang >= 14
# GCC >= 11
# CUDA >= 12.1 (Compute Capability 5.0)
# HIP >= 5.2
```

---

## 14.7 SIMD Auto-Detection

**File:** `cmake/gmxDetectSimd.cmake`

When `GMX_SIMD=AUTO`, CMake probes the CPU at configure time:

```cmake
function(gmx_suggest_simd _suggested_simd)
    gmx_run_cpu_detection(features)  # Run detection binary

    if(GMX_TARGET_X86)
        if(CPU_DETECTION_FEATURES MATCHES " avx512f ")
            # Intel: check FMA unit count (2 units → AVX-512, 1 → AVX2-256)
            gmx_detect_avx_512_fma_units(NUMBER_OF_AVX_512_FMA_UNITS)
            if(NUMBER_OF_AVX_512_FMA_UNITS EQUAL 2)
                set(OUTPUT_SIMD "AVX_512")
            else()
                set(OUTPUT_SIMD "AVX2_256")  # AVX2 faster with 1 FMA unit
            endif()
        elseif(CPU_DETECTION_FEATURES MATCHES " avx2 ")
            if(CPU_DETECTION_BRAND MATCHES "AMD")
                # Zen1: prefer 128-bit AVX2
                if("${CPU_DETECTION_FAMILY}" STREQUAL "23")
                    set(OUTPUT_SIMD "AVX2_128")
                else()
                    set(OUTPUT_SIMD "AVX2_256")
                endif()
            else()
                set(OUTPUT_SIMD "AVX2_256")
            endif()
        elseif(CPU_DETECTION_FEATURES MATCHES " avx ")
            if(CPU_DETECTION_FEATURES MATCHES " fma4 ")
                set(OUTPUT_SIMD "AVX_128_FMA")  # AMD Piledriver
            else()
                set(OUTPUT_SIMD "AVX_256")
            endif()
        elseif(CPU_DETECTION_FEATURES MATCHES " sse4.1 ")
            set(OUTPUT_SIMD "SSE4.1")
        elseif(CPU_DETECTION_FEATURES MATCHES " sse2 ")
            set(OUTPUT_SIMD "SSE2")
        endif()
    else()
        # ARM and IBM
        if(CPU_DETECTION_FEATURES MATCHES " sve ")
            set(OUTPUT_SIMD "ARM_SVE")
        elseif(CPU_DETECTION_FEATURES MATCHES " neon_asimd ")
            set(OUTPUT_SIMD "ARM_NEON_ASIMD")
        elseif(CPU_DETECTION_FEATURES MATCHES " vsx ")
            set(OUTPUT_SIMD "IBM_VSX")
        endif()
    endif()
endfunction()
```

### Smart Decisions

Notice the nuanced detection:
- **Intel AVX-512 with 1 FMA unit** (some Skylake): Uses AVX2-256 instead (higher throughput)
- **AMD Zen1**: Uses AVX2-128 (matches the internal 128-bit execution units)
- **AMD Piledriver**: Uses AVX-128-FMA (has FMA4 but no AVX2)

---

## 14.8 GPU Configuration

**Files:** `cmake/gmxManageCuda.cmake`, `cmake/gmxManageSycl.cmake`, `cmake/gmxManageOpenCL.cmake`, `cmake/gmxManageHip.cmake`

### GPU FFT Library Selection

```cmake
# CMakeLists.txt (lines 330-369)
# Maps GPU backend to default FFT library:
#   CUDA  → cuFFT
#   HIP   → VkFFT (or rocFFT)
#   OpenCL → VkFFT (or clFFT)
#   SYCL (ACPP) → VkFFT
#   SYCL (DPCPP) → MKL
```

### CUDA Architecture Selection

```cmake
# CUDA Compute Capability targets
# Default: auto-detect from installed GPU
# Can be overridden with: -DGMX_CUDA_ARCHITECTURES="70;80;90"
set_target_properties(target PROPERTIES CUDA_ARCHITECTURES "${GMX_CUDA_ARCHITECTURES}")
```

---

## 14.9 MPI Configuration

**File:** `cmake/gmxManageMPI.cmake`

```cmake
if (GMX_MPI)
    # Disable thread-MPI (incompatible)
    set(GMX_THREAD_MPI OFF CACHE BOOL "..." FORCE)
    set(GMX_LIB_MPI 1)

    find_package(MPI COMPONENTS CXX)
    if (NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "MPI support requested but no suitable compiler found")
    elseif (MPI_CXX_VERSION VERSION_LESS 3.0)
        message(FATAL_ERROR "MPI version 3.0 or higher required")
    endif()
else()
    set(GMX_LIB_MPI 0)
endif()
```

### Thread-MPI vs Library MPI

| Feature | Thread-MPI | Library MPI |
|---------|-----------|-------------|
| CMake option | `GMX_THREAD_MPI=ON` (default) | `GMX_MPI=ON` |
| Implementation | Threads within one process | Separate MPI processes |
| Multi-node | No | Yes |
| GPU-direct | No | Yes (with CUDA-aware MPI) |
| Typical use | Single workstation | HPC cluster |

---

## 14.10 FFT Library Configuration

**File:** `cmake/gmxManageFFTLibraries.cmake`

```cmake
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    if(GMX_DOUBLE)
        set(FFTW "FFTW")     # libfftw3 (double)
    else()
        set(FFTW "FFTWF")    # libfftw3f (single)
    endif()

    if(GMX_BUILD_OWN_FFTW)
        add_subdirectory(src/external/build-fftw)
    else()
        find_package(FFTW COMPONENTS ${LOWERFFTW})
        # Warn if FFTW lacks SIMD support
        if ((GMX_SIMD_ACTIVE MATCHES "SSE|AVX") AND NOT ${FFTW}_HAVE_SIMD)
            message(WARNING "FFTW without SIMD support — PME will be slow!")
        endif()
    endif()
elseif(${GMX_FFT_LIBRARY} STREQUAL "MKL")
    find_path(MKLROOT "include/mkl.h" PATHS ENV MKLROOT)
    # Link mkl_core, mkl_sequential, mkl_intel_lp64
endif()
```

**Important:** FFTW should be compiled with the same SIMD support as GROMACS. A FFTW without AVX on an AVX system will bottleneck PME.

---

## 14.11 Test Infrastructure CMake

**File:** `src/testutils/TestMacros.cmake`

```cmake
function (gmx_add_unit_test NAME EXENAME)
    cmake_parse_arguments(ARG "HARDWARE_DETECTION" "" "" ${ARGN})
    gmx_add_gtest_executable(${EXENAME} ${ARGN})
    gmx_register_gtest_test(${NAME} ${EXENAME} ${ARGN})
endfunction()

function (gmx_add_gtest_executable EXENAME)
    # Handles per-backend GPU compilation:
    if (GMX_GPU_CUDA)
        set_source_files_properties(${ARG_GPU_CPP_SOURCE_FILES}
            PROPERTIES LANGUAGE CUDA)
    elseif (GMX_GPU_HIP)
        set_source_files_properties(${ARG_GPU_CPP_SOURCE_FILES}
            PROPERTIES LANGUAGE HIP)
    elseif (GMX_GPU_SYCL)
        add_sycl_to_target(TARGET ${EXENAME}
            SOURCES ${ARG_SYCL_CPP_SOURCE_FILES} ${ARG_GPU_CPP_SOURCE_FILES})
    endif()

    target_link_libraries(${EXENAME} PRIVATE testutils gmock libgromacs)
endfunction()
```

---

## 14.12 SIMD Backend Comparison Table

| Backend | Width (float) | FMA | rsqrt bits | Platforms |
|---------|:---:|:---:|:---:|-----------|
| Reference | configurable | Emulated | Full | All (testing) |
| SSE2 | 4 | Emulated | 12 | x86 (2001+) |
| SSE4.1 | 4 | Emulated | 12 | x86 (2006+) |
| AVX-128-FMA | 4 | Hardware | 12 | AMD Piledriver |
| AVX-256 | 8 | Emulated | 12 | Intel Sandy Bridge |
| AVX2-128 | 4 | Hardware | 12 | AMD Zen1 |
| AVX2-256 | 8 | Hardware | 12 | Intel Haswell+ |
| AVX-512 | 16 | Hardware | 14 | Intel Skylake-X+ |
| ARM NEON | 4 | Hardware | — | ARM v8+ |
| ARM SVE | variable | Hardware | — | ARM v9 (Fugaku) |
| IBM VSX | 4 | Hardware | — | POWER8+ |

---

## 14.13 Exercises

### Exercise 1: Compare SIMD Backends
Open these three files:
- `impl_reference/impl_reference_simd_float.h`
- `impl_x86_sse2/impl_x86_sse2_simd_float.h`
- `impl_x86_avx2_256/impl_x86_avx2_256_simd_float.h`

Compare the `fma()` implementation. Count the instructions used in each.

### Exercise 2: Trace SIMD Detection
Run CMake on your system:
```bash
cmake -DGMX_SIMD=AUTO ..
```
Find the "Detected best SIMD instructions" message. Then:
1. Verify by checking your CPU's SIMD features (`lscpu` or `/proc/cpuinfo`)
2. Would a different SIMD setting be better for your CPU?

### Exercise 3: Newton-Raphson Precision
In `simd_math.h`, the `invsqrt()` function uses conditional Newton-Raphson iterations:
1. For SSE2 (`GMX_SIMD_RSQRT_BITS=12`), how many iterations are needed for single precision (24 bits)?
2. For AVX-512 (`GMX_SIMD_RSQRT_BITS=14`), how many iterations are needed?
3. Why does GROMACS check `RSQRT_BITS * 2 < ACCURACY_BITS` rather than `RSQRT_BITS < ACCURACY_BITS`?

### Exercise 4: Build Configuration
List all the CMake options you would use to build GROMACS for:
1. A laptop with Intel Core i7-12700H (AVX2)
2. An AMD EPYC 7763 cluster node with NVIDIA A100 GPU
3. A Fujitsu A64FX (ARM SVE) node

### Exercise 5: SIMD Kernel Selection
In `src/gromacs/nbnxm/kerneldispatch.cpp`:
1. Find where the kernel type is determined (4x4, 4xM, 2xMM)
2. What determines the choice between 4xM and 2xMM?
3. How does `GMX_SIMD_FLOAT_WIDTH` affect the cluster size M?

---

## 14.14 Key Takeaways

1. **Portable SIMD abstraction** — same code works across 10+ instruction sets via compile-time dispatch
2. **SimdFloat wraps hardware types** — `__m128` (SSE), `__m256` (AVX), `__m512` (AVX-512), `float32x4_t` (NEON)
3. **FMA is the key differentiator** — hardware FMA gives both speed and precision benefits
4. **Newton-Raphson iterations** refine hardware `rsqrt` estimates to full precision
5. **SIMD math uses minimax polynomials** — `log2`, `exp`, `erfc` are approximated with optimized coefficients
6. **CMake auto-detects optimal SIMD** — including nuances like Intel 1-FMA vs 2-FMA AVX-512 and AMD Zen1 128-bit preference
7. **Four GPU backends** configured through `GMX_GPU` option with backend-specific FFT library selection
8. **Thread-MPI vs Library MPI** — seamless switch between single-node and cluster modes
9. **FFTW SIMD matters** — FFTW should match GROMACS's SIMD level for PME performance
