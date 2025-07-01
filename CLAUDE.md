# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is AMReX

AMReX is a software framework for building massively parallel block-structured adaptive mesh refinement (AMR) applications. It provides infrastructure for solving partial differential equations on hierarchical adaptive grid structures with support for particles, embedded boundaries, and multiple parallelization backends (MPI, OpenMP, CUDA, HIP, SYCL).

## Build Commands

AMReX supports three build systems. The most commonly used are GNU Make and CMake.

### GNU Make Build
```bash
# Basic build with parallelism
make -j

# Common build configurations
make -j USE_MPI=TRUE USE_OMP=TRUE    # MPI + OpenMP build
make -j USE_CUDA=TRUE                 # GPU build with CUDA
make -j DEBUG=TRUE                    # Debug build
make -j DIM=2                         # 2D build (default is 3D)

# Clean commands
make clean         # Clean current build
make realclean     # Deep clean
```

### CMake Build
```bash
# Configure and build
cmake -S . -B build -DAMReX_MPI=YES -DAMReX_OMP=YES
cmake --build build -j

# Common CMake options
-DAMReX_SPACEDIM=3           # Set dimensionality (1, 2, or 3)
-DAMReX_GPU_BACKEND=CUDA     # Enable CUDA
-DAMReX_EB=YES               # Enable Embedded Boundary
-DAMReX_PARTICLES=YES        # Enable particles
```

## Running Tests

### Individual Test Example
```bash
cd Tests/LinearSolvers/ABecLaplacian_C/
make -j
./main3d.gnu.MPI.ex inputs
```

### Regression Testing
```bash
# Run full regression suite (requires regression_testing repo)
python regtest.py AMReX-tests.ini

# Run single test
python regtest.py --single_test MLMG_PoisLev AMReX-tests.ini
```

## Code Architecture

AMReX is organized into several core components under `Src/`:

- **Base/**: Fundamental data structures (`Box`, `MultiFab`, `ParallelDescriptor`)
- **AmrCore/**: Core AMR functionality (`AmrMesh`, `FillPatchUtil`, `FluxRegister`)
- **Amr/**: Full AMR framework with time stepping (`AmrLevel`, `StateData`)
- **LinearSolvers/MLMG/**: Multi-level multigrid solvers
- **Particle/**: Particle data structures and algorithms
- **EB/**: Embedded boundary (cut cell) support
- **Extern/**: Interfaces to external libraries (HYPRE, PETSc, SUNDIALS)

Key classes to understand:
- `Box`: Rectangular region in index space
- `BoxArray`: Collection of boxes (domain decomposition)
- `MultiFab`: Distributed array of `FArrayBox` objects (main data container)
- `Geometry`: Physical domain and coordinate system information
- `AmrMesh`: Manages hierarchical grid structure
- `ParticleContainer`: Manages particle data

## Development Practices

- Development happens on the `development` branch
- Pull requests should target `development`, not `main`
- Code style: 4 spaces indentation, space after function names in declarations
- Member variables prefixed with `m_`
- Use `amrex.abort_on_unused_inputs=1` when testing to catch unused parameters

## Common Development Tasks

### Adding a New Test
1. Create directory under `Tests/`
2. Add `GNUmakefile` that includes `$(AMREX_HOME)/Tools/GNUMake/Make.defs`
3. Write test in `main.cpp`
4. Add to regression test suite in `AMReX-tests.ini`

### Debugging
- Use `DEBUG=TRUE` in GNU Make builds
- Set `amrex.v=1` for verbose output
- Use `amrex.fpe_trap_invalid=1` to trap floating point errors
- `MultiFab::contains_nan()` to check for NaNs

### GPU Development
- Use `amrex::ParallelFor` for GPU kernels
- Mark device functions with `AMREX_GPU_DEVICE`
- Use `Gpu::DeviceVector` instead of `std::vector` in GPU code