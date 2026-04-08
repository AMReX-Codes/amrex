# Pipelined 2D-Slice TMA Backend for `StencilFor`

## Overview

This note describes a possible phase-2 CUDA backend for
`amrex::StencilFor` that uses true Tensor Memory Accelerator (TMA)
transfers on Hopper-class and newer NVIDIA GPUs.

The recommendation is deliberately narrow:

- do not replace the current row-wise `cub::BlockLoadToShared` path
- do not load the full 3-D staged cube in one bulk transfer
- do add a CUDA-specialized backend that pipelines 2-D slice transfers in
  the slow `z` direction

This follows the guidance from NVIDIA's GTC 2026 memory-bandwidth talk:
their 3-D stencil example reported a regression when loading the whole
3-D cube at once, and the improvement came from pipelining 2-D slices
through shared memory instead.

## Why A Separate Design Is Needed

The current CUDA implementation in `Src/Base/AMReX_StencilForG.H` stages
the full grown tile into shared memory and only starts compute after the
copy is complete.

That structure fits:

- the generic manual cooperative load path
- the current row-wise `BlockLoadToShared` path

It does not match the execution model that made TMA effective in the
stencil case shown in the NVIDIA slides:

- keep multiple slices in flight
- overlap later loads with current compute
- march through the tile in `z`
- recycle shared-memory slices as soon as the corresponding compute is done

For that reason, a true TMA backend should be treated as a new internal
execution strategy, not as a drop-in replacement for the existing
row-copy staging loop.

## Recommended Scope

The first pipelined TMA backend should stay narrower than the generic
`StencilFor` abstraction and fall back aggressively when the requested
case falls outside the tuned path.

Recommended initial scope:

- CUDA only
- `sm_90` and newer
- one staged component per launch
- cell-centered or otherwise regular packed `Array4` source layout
- uniform tiles in all dimensions
- full grown tile available in the source box
- `halo.z == 1` for the first implementation
- `TensorCopyPolicy::Auto` or `Always`

Recommended initial exclusions:

- partial tiles
- `halo.z > 1`
- multi-component staging
- user-controlled shared-padding modes that change the slice layout
- HIP, SYCL, or CPU specialization

The generic fallback and current row-wise async-copy path should remain
the default for all excluded cases.

## Why Not Full-Cube 3-D TMA

The full-cube approach has several drawbacks for this API:

- it preserves the current load-then-compute schedule, so it gives up the
  overlap that makes TMA attractive
- it increases the shared-memory footprint to the whole grown tile
- it makes occupancy and launch-shape tuning harder
- it does not match the best-performing pattern shown in the stencil
  slides

Inference from the NVIDIA material: if AMReX adds true TMA here, the
first target should be a pipelined slice backend, not a monolithic 3-D
cube transfer.

## Proposed Execution Model

The block still owns one logical 3-D compute tile, but it no longer
stages the full tile at once.

Instead:

1. Build one interior compute tile as today.
2. Grow it only in `x` and `y` for each staged slice.
3. Keep a ring buffer of shared-memory slices.
4. Prime the ring buffer with enough `z` planes to cover the current
   stencil window plus a small prefetch distance.
5. March through the tile in `z`, computing one or a few planes at a
   time.
6. As each oldest slice becomes dead, recycle that slot and issue the
   next TMA slice copy.

For a radius-1 stencil, a reasonable first shape is:

- live stencil window: `3` slices
- extra prefetch depth: `1` or `2` slices
- total ring size: `4` or `5` slices

This is much smaller than staging the full grown 3-D tile and creates a
natural place to overlap transfer and compute.

## Shared-Memory Layout

Each stage in the ring buffer is one 2-D `x-y` slice including `x` and
`y` halos.

Recommended layout rules:

- align each stage base to at least `128` bytes for ND TMA
- keep one shared-memory barrier per live stage slot
- reserve any bank-conflict padding as an internal backend detail
- avoid exposing stage padding or swizzle controls in the public API

The current `TileView` can still be used to expose the staged data to the
user lambda, but the backend will need an internal helper that remaps the
logical `z` coordinate onto the ring-buffer slot.

## Descriptor and Launch Model

This backend would need true CUDA TMA descriptors rather than the current
row-copy helper abstraction.

Recommended host-side responsibilities:

- cache `CUtensorMap` descriptors for the source layout
- key the cache on source pointer identity, source extents, source
  strides, datatype, tile shape, and halo shape
- pass descriptors to the kernel as immutable launch data

Recommended kernel-side responsibilities:

- elect one issuing thread or one issuing warp
- initialize the shared-memory barriers for the live stage slots
- issue TMA loads for the initial slices
- wait only on the stage needed for the next compute step
- recycle stage slots after all threads are done consuming them

The first implementation should use a single elected issuing thread.
Warp-specialized issuing can be revisited only if profiling shows that
instruction issue for TMA becomes a measurable bottleneck.

## Gating and Fallback Rules

The pipelined TMA backend should only activate when all of the following
are true:

- device capability is `sm_90+`
- tile shape is uniform
- `MT` is valid for the chosen compute mapping
- source base pointer is TMA-aligned
- source stride in bytes satisfies the TMA alignment rules
- staged slice byte count is a multiple of `16`
- shared-memory stage bases satisfy the stronger TMA alignment
- the ring-buffer footprint fits within the shared-memory budget
- the source box fully contains the grown tile region

If any one of these checks fails, the implementation should fall back in
this order:

1. row-wise `BlockLoadToShared` path when supported
2. manual cooperative shared-memory load
3. direct alias fallback

`TensorCopyPolicy::Always` should only mean "must use the selected native
backend". If the pipelined TMA path is not eligible, the call should
assert rather than silently dropping to another strategy.

## API Impact

The recommended design does not require a public API change.

The existing `StencilInfo` and `StencilFor` entry points can stay
unchanged if the TMA backend is treated as an internal specialization
selected from `TensorCopyPolicy::Auto`.

Possible internal additions:

- a CUDA-only pipeline state helper
- a CUDA-only slice descriptor cache
- a CUDA-only staged-tile accessor that maps logical `z` into ring slots

No public descriptor types should leak into the user-facing API.

## Testing Plan

Correctness coverage should compare the pipelined TMA backend against the
existing baseline stencil path and the current row-wise async-copy path.

Minimum correctness cases:

- 3-D 7-point Laplacian
- 3-D 27-point stencil with `halo = IntVect(1)`
- at least one tile shape with small `tile.z`
- at least one tile shape with larger `tile.z`
- forced fallback on unsupported shapes
- `TensorCopyPolicy::Always` failure on unsupported shapes

Performance validation should measure:

- runtime
- achieved bandwidth
- shared-memory footprint
- active blocks per SM
- sensitivity to pipeline depth
- sensitivity to `tile.z`

The comparison target should be the current row-wise async-copy backend,
not only the manual fallback.

## Open Questions

1. Should the first TMA slice backend support only `halo.z == 1`, or is
   there a small-extension path for `halo.z == 2` that is still worth the
   added shared-memory pressure?
2. Is one issuing thread sufficient, or does one issuing warp matter on
   Blackwell?
3. Do we need a swizzled shared-memory layout to avoid bank conflicts for
   wider `x-y` slices?
4. Where should descriptor caches live so their lifetime and invalidation
   rules remain obvious?
5. Should the backend require `require_full_tile = true` in the first
   implementation, even if other `StencilFor` backends allow partial-tile
   fallback?

## Recommendation

AMReX should not switch the current `StencilFor` CUDA path directly to
true 3-D TMA copies.

If AMReX pursues a second CUDA tensor-copy backend, it should be a
pipelined 2-D-slice design that:

- keeps the current CUB row-copy backend as the broadly applicable path
- targets only the cases where slice-pipelined TMA is likely to win
- uses strict eligibility checks
- proves value with dedicated Hopper and Blackwell benchmarks before the
  scope is widened

## References

- NVIDIA CUDA Programming Guide, asynchronous data copies:
  <https://docs.nvidia.com/cuda/cuda-programming-guide/04-special-topics/async-copies.html>
- NVIDIA GTC 2026 slide deck:
  <https://github.com/user-attachments/files/26574763/S81666.pdf>
- Existing AMReX stencil prototype note:
  `Docs/Notes/TensorCopyStencilPrototype.md`
