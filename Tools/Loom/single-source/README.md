# Single-source LOOM build of AMReX

This directory reproduces the single-source LOOM test: AMReX's single-GPU tests
built with one C++ frontend (cxx) for both host and device code, and one set of
libc++ headers for both. cxx compiles each file twice:

1. to x86-64 host code, through cxx's MLIR/LLVM back end;
2. through Loom's C++ importer (which embeds the same cxx) and `loom-compile`,
   to a gfx90a code object that is embedded in the host object.

HRX is the GPU runtime, on its native AMD GPU driver (libamdf); no part of
ROCm is used. The pass criterion is that `GPU_ArenaLifetime_3d`,
`GPU_CNS_Sod_3d` and `GPU_ExternalStream_3d` pass, and that the CNS run agrees
with `cns-expected.txt` to round-off.

The device pass still sees a prelude (`../amrex_loom_device_prelude.H`)
and Loom's predefined `__loom__`, which the host pass does not, so the build is
not yet strictly single-source. Only the backend header
`../hip/hip/hip_runtime.h` tests `__loom__`, to define
`AMREX_LOOM_DEVICE_PASS` for the rest of the backend; no AMReX source under
`Src/` tests either macro.

## Requirements

- Linux x86-64 with an AMD gfx90a GPU (MI210, MI250X).
- The amdgpu kernel driver with access to `/dev/kfd`, which HRX's libamdf uses
  directly. ROCm is not needed: no ROCr/HSA runtime, HIP, device libraries or
  ROCm headers.
- GCC 14 (its install directory, for libstdc++, glibc crt files and compiler
  headers).
- glibc headers (`/usr/include` on most systems).
- [pixi](https://pixi.sh) for the conda-forge packages in `pixi.toml`:
  clang 22.1.8 (builds the tools and links), MLIR/LLVM 22.1.8 (cxx's back end),
  libc++ 23.1.5 (host and device), CMake and Ninja.
- About 30 GB of disk, and a few GB of memory per build job.

## Steps

```bash
mkdir loom-work && cd loom-work
git clone -b loom-single-source https://github.com/BenWibking/cplusplus.git cplusplus
git clone -b loom-single-source https://github.com/BenWibking/hrx.git hrx
git clone -b loom-single-source https://github.com/BenWibking/amrex.git amrex
git -C cplusplus checkout --detach eb342fe42f   # pinned cxx commit
git -C hrx checkout --detach 2ef3663193         # pinned HRX commit

S=$PWD/amrex/Tools/Loom/single-source
pixi install --manifest-path $S/pixi.toml   # conda-forge toolchain pieces
```

Set the machine-specific paths at the top of `env.sh`, or export them before
running the scripts:

| Variable          | Meaning                                  | Default                                  |
|-------------------|------------------------------------------|------------------------------------------|
| `GCC_INSTALL_DIR` | GCC 14 `lib/gcc/<triple>/<version>`      | `/usr/lib64/gcc/x86_64-suse-linux/14`    |
| `GLIBC_INCLUDE`   | glibc headers                            | `/usr/include`                           |
| `AMREX_LOOM_SYSROOT` | sysroot for host compiles and links; same glibc as `GLIBC_INCLUDE` | `/` |
| `LOOM_WORK`       | directory holding the three clones       | current directory                        |
| `BUILD_JOBS`      | parallel build jobs                      | `8`                                      |

Then, from `loom-work` (on a GPU node for the last step):

```bash
$S/build-cxx.sh     # cxx with its MLIR back end              (build/cxx)
$S/build-hrx.sh     # loom-import-cxx, loom-compile, libhrx    (build/hrx)
$S/build-amrex.sh   # AMReX LOOM backend and the three tests   (build/amrex)
$S/run-tests.sh     # ctest for the three tests, then the CNS check
```

`build-hrx.sh` points HRX's cxx dependency at the `cplusplus` clone
(`FETCHCONTENT_SOURCE_DIR_CXX`), so Loom's importer and the host compiler are
built from the same cxx source. It builds only HRX's native AMD GPU driver, and
`env.sh` selects it (`HRX_GPU_DRIVER=amdgpu`).

## Expected result

```text
100% tests passed, 0 tests failed out of 3
steps checked: 20; max relative difference from the recorded run: ... (tolerance 1e-10); bit-identical: ...
  mass: 147456.0 (CPU reference 147455.99999999956, relative 3.0e-15)
  ...
PASS
```

On a different GPU model the CNS numbers may differ from the recorded MI210 run
in the last digits, because GPU reductions can sum in a different order; the
check allows a relative difference of 1e-10 (`CNS_RTOL` overrides it).

## Notes

- Memory: keep `BUILD_JOBS` modest. Each cxx host compile and each device import
  takes a few GB.
- Executables link with `-no-pie`: cxx does not generate position-independent
  code yet.
- Kernels that the importer cannot translate are dropped and reported by
  `amrex-loom-cxx` (`dropped N launched kernel(s)`); none of them is reached by
  the three tests.
- The driver and the host-compiler wrapper are `../amrex-loom-cxx` and
  `../amrex-loom-host-cxx`; both document their environment variables.
