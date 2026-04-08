# Tensor Copy Stencil Prototype Plan

## Overview

This document describes a prototype implementation plan for adding a
C++ lambda abstraction to AMReX for block-tiled stencil kernels that can
stage multidimensional tiles from global memory into shared memory or LDS
using next-generation tensor copy instructions where available.

The motivating discussion is AMReX issue
[`#5279`](https://github.com/AMReX-Codes/amrex/issues/5279), which calls
out:

- NVIDIA Tensor Memory Accelerator (TMA) on Hopper (`sm_90`) and newer
- AMD Tensor Data Mover (TDM) on `gfx1250` and newer
- the need for an AMReX abstraction suitable for multidimensional stencil
  kernels

The proposed prototype is intentionally narrow:

- add one new opt-in block-stencil launch path
- keep existing `ParallelFor` and `launch` unchanged
- provide a generic fallback path on all backends
- provide a first native tensor-copy implementation on CUDA/Hopper-class
  hardware
- defer the first native AMD implementation until the API shape is proven
- explicitly support wide-halo stencils with radius greater than one,
  including practical cases with halo widths of 4 or 5 cells

## Motivation

AMReX already has efficient pointwise kernel abstractions built around
`ParallelFor`, and more general block launches via `amrex::launch`.
Those are a good fit for many kernels, but they do not directly capture
the execution model of tensor copy instructions:

- one thread block owns one compute tile
- the block cooperatively stages a grown tile into shared memory or LDS
- the block waits for the staged copy to complete
- the block computes stencil values over the interior tile

That execution model is fundamentally different from AMReX's default
flattened pointwise `ParallelFor` path.

The prototype should therefore introduce a dedicated block-stencil
abstraction instead of overloading pointwise `ParallelFor` further.

## Current AMReX Constraints

This plan is based on the current AMReX implementation.

### 1. `ParallelFor` is pointwise

`ParallelFor` dispatches pointwise lambdas over a flattened cell space.
The lambda may optionally accept `Gpu::Handler`, but the launch model is
still one logical cell per iteration.

Relevant files:

- `Src/Base/AMReX_GpuLaunchFunctsG.H`
- `Src/Base/AMReX_GpuLaunchFunctsC.H`
- `Docs/sphinx_documentation/source/GPU.rst`

### 2. `KernelInfo` is too small for this feature

`Gpu::KernelInfo` currently only carries a reduction flag. It is not the
right place to encode tile shape, halo size, copy policy, padding, or
backend capability selection.

Relevant file:

- `Src/Base/AMReX_GpuKernelInfo.H`

### 3. `amrex::launch` is the better starting point

AMReX already uses `amrex::launch(..., shared_mem_bytes, ..., lambda)` for
block-cooperative shared-memory kernels. This is the natural substrate
for a new stencil abstraction.

Relevant files:

- `Src/Base/AMReX_BaseFabUtility.H`
- `Src/Base/AMReX_GpuContainers.H`
- `Docs/sphinx_documentation/source/GPU.rst`

### 4. `ExecutionConfig(Box)` assumes 1-D decomposition

`Gpu::ExecutionConfig(Box)` currently uses a 1-D decomposition of the cell
space, and the source comment notes that this assumption matters for
existing code. A tensor-copy stencil path should not reuse that mapping
for the compute phase. It should launch one block per tile instead.

Relevant file:

- `Src/Base/AMReX_GpuLaunch.H`

### 5. `Array4` is packed, not explicitly strided

`Array4` carries the source pointer, bounds, and packed stride
information. That is useful for source metadata, but the current
constructor path computes packed strides from the extents. It does not
directly represent a padded shared-memory tile layout.

That means the prototype should not assume that a padded tile can be
represented as a plain `Array4`.

Relevant file:

- `Src/Base/AMReX_Array4.H`

### 6. Backend specialization already exists elsewhere

AMReX already wraps backend-specific block collectives through
backend-specialized implementations using CUB, rocPRIM, and SYCL
facilities. The tensor-copy abstraction should follow the same pattern.

Relevant files:

- `Src/Base/AMReX_GpuReduce.H`
- `Src/Base/AMReX_Scan.H`

### 7. Alignment must be checked explicitly

AMReX arena allocation currently documents a default arena alignment of
16 bytes. Tensor-copy paths will have stricter backend- and datatype-
dependent alignment requirements, so the prototype must validate those
requirements and fall back automatically when they are not met. In
particular, CUDA async-copy helpers such as `cub::BlockLoadToShared`
require not only an aligned base shared-memory pointer but also aligned
per-copy destinations, so row strides derived from `shared_padding`
must preserve the required shared-memory buffer alignment.

Relevant file:

- `Src/Base/AMReX_Arena.H`

## Prototype Goals

The prototype should:

- preserve existing AMReX launch APIs
- introduce one new opt-in abstraction for stencil kernels
- hide all backend-specific descriptor and instruction details from user code
- support a generic fallback path on CPU, HIP, SYCL, and older CUDA GPUs
- support a native CUDA tensor-copy path on Hopper-class GPUs
- keep the user kernel in ordinary C++ lambda form
- work with existing `Array4` global-memory inputs
- support halo widths from 1 up to at least 5 cells in each dimension
- expose enough tile metadata to support performance tuning
- make fallback behavior explicit and predictable

## Non-Goals

This prototype does not try to:

- redesign `ParallelFor`
- generalize to every block-structured kernel shape immediately
- stage multiple source arrays in the first implementation
- support tensor-store instructions in the first implementation
- expose vendor descriptor types in the public API
- land a full AMD native tensor-copy path in the first patch series
- change the storage layout of `FArrayBox` or `MultiFab`

## Proposed Public API

The prototype should add a new launch family rather than adding more flags
to `ParallelFor`.

### High-level launch object

```cpp
namespace amrex::Gpu {

enum class TensorCopyPolicy {
    Auto,
    Always,
    Never
};

struct StencilInfo {
    IntVect tile = IntVect(AMREX_D_DECL(32, 8, 4));
    IntVect halo = IntVect(1);
    IntVect shared_padding = IntVect(0);
    TensorCopyPolicy tensor_copy = TensorCopyPolicy::Auto;
    bool require_full_tile = false;
    bool stage_one_component = true;

    StencilInfo& setTile (IntVect const& v) noexcept;
    StencilInfo& setHalo (IntVect const& v) noexcept;
    StencilInfo& setSharedPadding (IntVect const& v) noexcept;
    StencilInfo& setTensorCopyPolicy (TensorCopyPolicy v) noexcept;
    StencilInfo& setRequireFullTile (bool v) noexcept;
};

}
```

### High-level launch functions

Prototype entry points:

```cpp
template <int MT = AMREX_GPU_MAX_THREADS, typename T, typename F>
void StencilFor (Box const& bx,
                 Gpu::StencilInfo const& info,
                 Array4<T const> const& src,
                 F&& f) noexcept;

template <int MT = AMREX_GPU_MAX_THREADS, typename T, typename F>
void StencilFor (Box const& bx,
                 int ncomp,
                 Gpu::StencilInfo const& info,
                 Array4<T const> const& src,
                 F&& f) noexcept;
```

The staged source field is explicit. Destination arrays and any other
captured state stay in the user lambda.

Example target usage:

```cpp
auto info = amrex::Gpu::StencilInfo{}
    .setTile(IntVect(AMREX_D_DECL(32, 8, 4)))
    .setHalo(IntVect(1))
    .setTensorCopyPolicy(amrex::Gpu::TensorCopyPolicy::Auto);

amrex::StencilFor<128>(bx, info, phi_old,
[=] AMREX_GPU_DEVICE (auto const& tile, int i, int j, int k) noexcept
{
    phi_new(i,j,k) = c0 * tile(i,j,k)
                   + c1 * (tile(i-1,j,k) + tile(i+1,j,k)
                         + tile(i,j-1,k) + tile(i,j+1,k)
                         + tile(i,j,k-1) + tile(i,j,k+1));
});
```

This keeps the user-visible stencil expression close to existing
`Array4`-based kernels while giving AMReX ownership of the tile staging
step.

## Accessor Design

The staged tile accessor should be a new lightweight view type rather than
plain `Array4`.

Recommended prototype type:

```cpp
template <typename T>
class TileView {
public:
    AMREX_GPU_DEVICE T& operator() (int i, int j, int k) const noexcept;
    AMREX_GPU_DEVICE T& operator() (int i, int j, int k, int n) const noexcept;

    AMREX_GPU_DEVICE bool contains (int i, int j, int k) const noexcept;
    AMREX_GPU_DEVICE Box tileBox () const noexcept;
    AMREX_GPU_DEVICE Box validBox () const noexcept;
};
```

`halo` must not be treated as a special-case radius-1 parameter. The
prototype should support arbitrary positive halo widths, with validation
and testing covering at least the range `[1,5]`.

Design choices:

- index with absolute cell coordinates, not relative offsets
- carry explicit strides so padded shared-memory layouts are representable
- keep the staged region's grown box available for debug checks
- optionally expose component indexing for a later multi-component path

This avoids changing `Array4` for the first prototype while still allowing
padding and backend-specific shared-memory layout rules.

## Internal Execution Model

### 1. Tile decomposition on the host

For each `Box bx` passed by the user:

- decompose `bx` into compute tiles using `StencilInfo::tile`
- assign one thread block per compute tile
- define the staged tile as `compute_tile.grow(info.halo)`
- intersect that staged tile with the available source region
- if the staged tile cannot be satisfied and `require_full_tile == true`,
  abort in debug or fall back in optimized builds

This launch should not reuse `ParallelFor(Box)` because the block-to-tile
ownership is the key abstraction.

Wide halos directly affect tile feasibility. A tile size that is sensible
for `halo = 1` may become impossible for `halo = 4` or `5`, so the launch
path must validate the staged tile footprint before selecting the native
copy path.

### 2. Shared-memory layout construction

For each compute tile:

- compute interior extents = `tile`
- compute staged extents = `tile + 2 * halo`
- apply optional shared-memory padding
- allocate one dynamic shared-memory buffer per block
- construct `TileView` over that buffer with explicit strides

The first prototype should support only one staged component per launch.
That keeps descriptor construction, shared-memory sizing, and fallback
logic simple.

For wide halos, shared-memory use grows quickly. In 3-D, a radius-5 halo
adds 10 cells in every dimension, so staged extents can become much larger
than the interior tile. The prototype should therefore:

- compute shared-memory requirements on the host before launch
- reject native tensor-copy staging if the staged tile does not fit
- fall back automatically to a smaller-tile or generic path when needed

### 3. Copy phase

The block enters one of three paths:

1. native tensor-copy path
2. optimized cooperative block-load path
3. generic manual fallback path

#### Native path

Used only when:

- the backend reports tensor-copy support
- the device architecture is new enough
- the pointer, bounds, and strides satisfy alignment and size constraints
- each staged row start in shared memory remains aligned, not just the
  base pointer
- the staged tile shape fits within shared-memory limits
- the halo width is supported by the selected descriptor construction path

#### Cooperative load path

Used when the backend has a block-load helper but native tensor-copy is
not available or not selected.

#### Generic fallback path

Each thread block cooperatively loads the staged tile from global memory
using ordinary global loads and shared-memory stores.

This fallback path is especially important for wide-halo stencils because
some halo widths or tile shapes may be correct but not profitable, or may
exceed native tensor-copy constraints on a given backend.

### 4. Synchronization

After the load phase:

- wait for the copy operation to complete
- synchronize the block
- enter the compute phase

### 5. Compute phase

Threads iterate over the interior tile only. Each thread:

- maps its local thread id to one or more cells in the compute tile
- calls the user lambda with `TileView` and absolute cell indices
- writes results directly to globally captured destination arrays

The first prototype should not attempt a tensor-store path. A direct
global store is simpler and sufficient to validate the staging model.

## Backend Strategy

### CPU

The CPU path should preserve ordinary correctness and developer
ergonomics:

- allocate a small stack or heap tile buffer if convenient, or bypass tile
  staging entirely
- present the same lambda shape to the user
- prioritize clarity over aggressive optimization

The goal is API portability, not a CPU performance feature.

### CUDA

The first native implementation should be CUDA-first.

Why:

- TMA support is already documented in CCCL/libcudacxx
- the issue already references NVIDIA's `cub::BlockLoadToShared`
- AMReX already uses CUB in other backend-specialized code
- CUDA capability checks are already available through
  `Gpu::Device::devicePropMajor()` and `devicePropMinor()`

Recommended CUDA rollout:

#### Phase 1 CUDA native path

Use the highest-level supported CUDA abstraction that still maps to TMA.
Prefer:

- `cub::BlockLoadToShared` if it covers the needed layout cases

Fallback within CUDA:

- direct TMA descriptor construction using CCCL/libcudacxx if the higher-
  level helper does not fit AMReX's layout or padding needs

The prototype should be guarded on:

- `AMREX_USE_CUDA`
- toolkit and CCCL availability
- `Gpu::Device::devicePropMajor() >= 9`

#### Descriptor translation

If direct TMA descriptors are needed, AMReX should build them internally
from `Array4` metadata. The dimension ordering and stride conventions must
be translated carefully because AMReX uses x-fastest storage while the
CCCL TMA APIs are documented in DLTensor terms.

This translation must remain an internal implementation detail.

#### Descriptor caching

The prototype should start simple:

- create one descriptor per source array and launch configuration on the host
- pass it into the kernel launch

If the prototype is promising, a later optimization can cache descriptors
by a key such as:

- source pointer
- source extents and strides
- datatype
- component
- tile shape
- halo

Wide halos make descriptor caching more valuable because the descriptor
setup cost becomes a smaller fraction of total work as the staged region
gets larger.

### HIP / AMD

The first patch series should not block on a native AMD implementation.

Recommended HIP plan:

- Phase 1: generic fallback path only
- Phase 2: add `supportsTensorCopy()` detection for `gfx1250+`
- Phase 3: add a native TDM-backed implementation once the relevant LLVM or
  HIP APIs are stable enough for AMReX to wrap cleanly

The public API should be designed now so that the later AMD path can slot
in without changing user code.

### SYCL

The prototype should support SYCL only through the generic fallback path.

There is no reason to block the API on a SYCL-native tensor-copy story.

## Capability and Fallback Interface

AMReX should expose one internal capability layer:

```cpp
namespace amrex::Gpu::detail {

struct TensorCopyCaps {
    bool available = false;
    bool requires_alignment_check = true;
    bool supports_multidim = false;
    bool supports_shared_padding = false;
};

TensorCopyCaps queryTensorCopyCaps () noexcept;

bool canUseTensorCopy (void const* p,
                       Box const& src_box,
                       IntVect const& tile,
                       IntVect const& halo,
                       std::size_t element_size) noexcept;
}
```

And optionally one public coarse-grained query:

```cpp
namespace amrex::Gpu {
bool supportsTensorCopy () noexcept;
}
```

This keeps the public surface small while still giving the implementation
room to make backend-specific decisions.

## File Layout

The prototype should be introduced in a way that matches AMReX's existing
header organization.

Recommended files:

- `Src/Base/AMReX_StencilFor.H`
  User-facing declarations and small wrappers.
- `Src/Base/AMReX_StencilForC.H`
  CPU implementation.
- `Src/Base/AMReX_StencilForG.H`
  GPU implementation and backend dispatch.
- `Src/Base/AMReX_GpuTensorCopy.H`
  Internal capability checks, descriptor wrappers, and backend hooks.
- `Src/Base/AMReX_Gpu.H`
  Add the new include once the API is ready.
- `Src/Base/CMakeLists.txt`
  Add headers to the build.
- `Src/Base/Make.package`
  Add headers to the GNUmake package list.

If `TileView` proves broadly useful, it can later move into its own header.
For the prototype, it should stay scoped to the stencil-launch feature.

## Recommended Implementation Phases

### Phase 0: scaffolding

- add the new headers
- add `StencilInfo`
- add a CPU fallback path
- add a GPU fallback path with manual cooperative loads
- add one focused GPU test

Deliverable:

- new API compiles and produces correct answers everywhere

### Phase 1: CUDA native path

- add CUDA capability detection
- add alignment checks
- integrate `cub::BlockLoadToShared` or direct TMA descriptors
- keep the manual fallback path active for unsupported cases

Follow-on note:

- `Docs/Notes/StencilForPipelinedTMA.md` describes a possible later CUDA
  backend based on pipelined 2-D slice TMA rather than row-wise block
  loads or full-cube 3-D TMA.

Deliverable:

- CUDA/Hopper can use the native path automatically

### Phase 2: performance validation

- benchmark 7-point and 27-point stencil kernels
- tune default tile sizes
- measure shared-memory footprint and occupancy tradeoffs
- validate fallback behavior on misaligned and unsupported inputs

Deliverable:

- evidence that the abstraction is worthwhile

### Phase 3: API hardening

- decide whether `TileView` should remain local or become a more general
  explicit-stride view type
- decide whether the 4-D staged-component overload should be enabled
- decide whether descriptor caching is needed

Deliverable:

- stable prototype API ready for broader use

### Phase 4: AMD native path

- add architecture detection for `gfx1250+`
- wrap the best available TDM interface
- validate against the same tests and benchmarks

Deliverable:

- native AMD tensor-copy implementation behind the same API

## Testing Plan

The prototype should include both correctness and performance checks.

### Correctness tests

Add a focused test such as:

- `Tests/GPU/TensorCopyStencil/main.cpp`

Test matrix:

- 2-D and 3-D
- halo widths `1`, `2`, `4`, and `5`
- both isotropic halos like `IntVect(4)` and at least one anisotropic case
  such as `IntVect(AMREX_D_DECL(5,2,1))`
- single component staged source
- fallback path forced on
- tensor-copy path forced on when supported
- edge tiles and partial tiles
- ghost-cell-dependent stencil values

The test kernels should include both:

- radius-1 stencils, to compare against the current common case
- wider-radius stencils that consume halo values at distances 2 through 5

This avoids accidentally validating only the descriptor plumbing while
failing to exercise the wider-halo indexing paths.

The result should be compared against a baseline `ParallelFor` stencil
implementation similar to the heat equation example in
`Tests/HeatEquation/main.cpp`.

### Performance tests

Measure:

- kernel time
- achieved bandwidth
- occupancy or active blocks per SM where available
- sensitivity to tile shape and block size

Initial benchmark kernels:

- 7-point Laplacian
- 27-point stencil
- one wide-halo stencil in 2-D
- one wide-halo stencil in 3-D

Suggested wide-halo benchmark shapes:

- 2-D radius-4 or radius-5 box stencil
- 3-D radius-2 or radius-4 structured stencil, depending on shared-memory
  footprint

### Negative tests

Explicitly test fallback behavior for:

- unsupported architectures
- misaligned pointers
- staged tile too large for shared memory
- tile-plus-halo region not fully available
- halo widths that force fallback because the staged region is too large for
  the selected native copy path

## Risks and Open Questions

### 1. Shared-tile accessor shape

The biggest immediate design question is whether to:

- introduce a narrow `TileView` for this feature, or
- generalize `Array4` to support explicit strides

The prototype should choose `TileView` first to avoid broad churn.

### 2. CUDA helper coverage

It is not yet guaranteed that `cub::BlockLoadToShared` will cover all
AMReX layout and padding cases needed for a production path. If it does
not, the CUDA implementation may need to drop to a lower-level TMA
descriptor interface sooner.

### 3. Multi-component staging

Staging multiple components in one shot is desirable, but it increases:

- descriptor complexity
- shared-memory footprint
- padding complexity
- fallback-path complexity

The prototype should stage one component per launch first.

This is even more important when the halo radius is 4 or 5, because wide
halos already put strong pressure on shared-memory capacity.

### 4. Tile-size defaults for wide halos

Reasonable default interior tile sizes will likely depend on halo width.
For example, a tile shape that is appropriate for radius-1 stencils may be
too large once the staged region includes a radius-5 halo.

The prototype should therefore avoid assuming one universal default tile.
At minimum, the implementation should:

- validate the user-requested tile against the shared-memory budget
- document that wide-halo stencils may need smaller interior tiles
- consider choosing backend defaults as a function of halo width in a later
  tuning pass

### 4. Descriptor cache lifetime

If AMReX caches native descriptors, the invalidation rules must be clear.
Source pointer identity alone may not be enough if the logical extents or
component slice change.

### 5. AMD implementation surface

The public AMD documentation currently points more clearly to compiler or
IR-level constructs than to a stable user-level C++ runtime API. That is a
strong reason to keep the AMD-native path out of the first patch series.

## Recommended Initial Patch Series

1. Add `StencilInfo`, `TileView`, `StencilFor`, and the generic fallback.
2. Add one correctness test and one benchmark-style test.
3. Add CUDA native staging behind `TensorCopyPolicy::Auto`.
4. Tune tile defaults and document fallback behavior.
5. Revisit AMD-native support after the CUDA path is validated.

This keeps the first series small enough to review while still proving the
core idea.

## References

### AMReX issue

- AMReX issue `#5279`: <https://github.com/AMReX-Codes/amrex/issues/5279>

### AMReX source references

- `Src/Base/AMReX_GpuKernelInfo.H`
- `Src/Base/AMReX_GpuLaunch.H`
- `Src/Base/AMReX_GpuLaunchFunctsG.H`
- `Src/Base/AMReX_GpuLaunchFunctsC.H`
- `Src/Base/AMReX_BaseFabUtility.H`
- `Src/Base/AMReX_GpuContainers.H`
- `Src/Base/AMReX_Array4.H`
- `Src/Base/AMReX_GpuReduce.H`
- `Src/Base/AMReX_Scan.H`
- `Src/Base/AMReX_Arena.H`
- `Docs/sphinx_documentation/source/GPU.rst`
- `Tests/HeatEquation/main.cpp`

### NVIDIA references

- CCCL Tensor Memory Accelerator overview:
  <https://nvidia.github.io/cccl/libcudacxx/extended_api/tma.html>
- CCCL `make_tma_descriptor`:
  <https://nvidia.github.io/cccl/unstable/libcudacxx/extended_api/tma/make_tma_descriptor.html>
- CUB `BlockLoadToShared` reference from the issue discussion:
  <https://github.com/NVIDIA/cccl/blob/f296704dbf65dc72bf09cc3ab587a650b9e61292/cub/cub/block/block_load_to_shared.cuh#L53>

### AMD references

- MLIR AMDGPU dialect:
  <https://mlir.llvm.org/docs/Dialects/AMDGPU/>
