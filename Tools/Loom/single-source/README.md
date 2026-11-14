# Single-source LOOM build of AMReX

This directory reproduces the single-source LOOM test: AMReX's single-GPU tests
built with one C++ frontend (cxx) for both host and device code, and one set of
libc++ headers for both. Each file is compiled by Loom's offload driver,
`loom-cxx`, which parses and type-checks it once with cxx and gives that one
translation unit to both back ends:

1. host code, through cxx's MLIR/LLVM back end (first, so device import cannot
   change it);
2. the launched kernels, through Loom's C++ importer and the Loom compiler, to
   a gfx90a code object that is embedded in the same object file, with the
   kernel descriptors the host code references.

HRX is the GPU runtime, on its native AMD GPU driver (libamdf); no part of
ROCm is used. The pass criterion is that `GPU_ArenaLifetime_3d`,
`GPU_CNS_Sod_3d`, `GPU_ExternalStream_3d` and `GPU_LoomCodeObjects_3d` pass,
and that the CNS run agrees with `cns-expected.txt` to round-off.
`GPU_LoomCodeObjects_3d` checks that the runtime loads a translation unit's
code object only at the first launch of one of its kernels, and that kernel
arguments narrower than a kernel argument slot reach the device.

Host and device code come from the same parse, so they see the same headers,
including Loom's `loomcxx/` headers. Loom's importer lowers the C math functions
and their `__builtin_*` forms itself, so no device prelude is injected.

Inside functions that run on both sides, `AMREX_IF_ON_DEVICE` and
`AMREX_IF_ON_HOST` expand to ordinary `if` statements on
`loom::target::is_device()` (`<loomcxx/target.h>`). The host compiler sees a
constant `false` and emits only the host arm, even at `-O0`; the importer
folds the call to `true` and imports only the device arm. Both arms are parsed
(once).

Device code is written on Loom's own C++ API, loomcxx
(`Src/Base/AMReX_LoomDevice.H`): kernel topology, workgroup barriers and
`[[loom::workgroup]]` storage, atomics, fences, subgroup shuffles and
`loom::scalar` math. Loom operations are declarations, which the host never
emits because only kernel bodies and device arms reach them. Host code names
kernel templates only through `loom::kernel::kernel_ref`. Nothing in AMReX
tests which back end it is compiled for.

## Requirements

- Linux x86-64 with an AMD gfx90a GPU (MI210, MI250X).
- The amdgpu kernel driver with access to `/dev/kfd`, which HRX's libamdf uses
  directly. ROCm is not needed: no ROCr/HSA runtime, HIP, device libraries or
  ROCm headers.
- GCC 14 (its install directory, for libstdc++, glibc crt files and compiler
  headers).
- glibc headers (`/usr/include` on most systems).
- [pixi](https://pixi.sh) for the conda-forge packages in `pixi.toml`:
  clang 22.1.8 (builds the tools and links), MLIR/LLVM 22.1.8 (cxx's back end,
  linked into `loom-cxx`), libc++ 23.1.5 (host and device), CMake and Ninja.
- About 30 GB of disk, and a few GB of memory per build job.

## Steps

```bash
mkdir loom-work && cd loom-work
git clone -b loom-single-source https://github.com/BenWibking/cplusplus.git cplusplus
git clone -b loom-single-source https://github.com/BenWibking/hrx.git hrx
git clone -b loom-single-source https://github.com/BenWibking/amrex.git amrex
git -C cplusplus checkout --detach 0caf64b6ee   # pinned cxx commit
git -C hrx checkout --detach a74773957c         # pinned HRX commit

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
$S/build-hrx.sh     # loom-cxx (with cxx and its MLIR back end), libhrx  (build/hrx)
$S/build-amrex.sh   # AMReX LOOM backend and the four tests              (build/amrex)
$S/run-tests.sh     # ctest for the four tests, then the CNS check
```

`build-hrx.sh` configures HRX with `LOOM_CXX_HOST_CODEGEN=ON` and points HRX's
cxx dependency at the `cplusplus` clone (`FETCHCONTENT_SOURCE_DIR_CXX`): HRX
builds cxx through cxx's own CMake project with its MLIR back end (pixi's
MLIR/LLVM), and `loom-cxx` links that one copy of cxx with Loom's importer and
compiler. pixi's MLIR package is used as installed. HRX builds only its native
AMD GPU driver, and `env.sh` selects it (`HRX_GPU_DRIVER=amdgpu`).

## Expected result

```text
100% tests passed, 0 tests failed out of 4
steps checked: 20; max relative difference from the recorded run: ... (tolerance 1e-10); bit-identical: ...
  mass: 147456.0 (CPU reference 147455.99999999956, relative 3.0e-15)
  ...
PASS
```

On a different GPU model the CNS numbers may differ from the recorded MI210 run
in the last digits, because GPU reductions can sum in a different order; the
check allows a relative difference of 1e-10 (`CNS_RTOL` overrides it).

## Notes

- `build-hrx.sh` builds only `loom-cxx` and libhrx. Loom's standalone
  importer and compiler, `loom-import-cxx` and `loom-compile`, are debugging
  tools: `cmake --build build/hrx --target loom-import-cxx loom-compile`.
- Memory: keep `BUILD_JOBS` modest. Each `loom-cxx` compile takes a few GB.
- Executables link with `-no-pie`: cxx does not generate position-independent
  code yet.
- Kernels that the importer cannot translate are dropped and reported by
  `loom-cxx` (`dropped N launched kernel(s)`, `--loom-drop-failed-kernels`);
  none of them is reached by the four tests. Host errors, and any other
  failure of the device pass, are compile errors.
- AMReX is configured with the toolchain file `../LoomToolchain.cmake`:
  `loom-cxx` is the C++ compiler, called once per compile, and `clang++` with
  libc++ links; the file documents its settings, which `env.sh` exports.
  AMReX's CMake adds the device options (`-DAMReX_AMD_ARCH` and the drop
  policy) and libhrx; `loom-cxx` finds
  Loom's `loomcxx/` headers next to itself, as clang finds its resource
  directory. `loom-cxx --loom-stats` reports, per file, that the source was
  preprocessed and parsed once, and the time of each phase.
