.. _sec:nondeterminism:

Floating-Point Nondeterminism
=============================

AMReX executes most kernels in parallel on CPUs, GPUs, and across MPI ranks.
Because floating-point addition and multiplication are not associative, the
final rounding depends on the precise order in which partial updates are
executed.  This section summarizes the places in the AMReX source that rely on
atomics, reductions, or collective communication whose scheduling can vary from
run to run, leading to small changes in floating-point answers even when inputs
match.

Atomic accumulations in FAB-based operations
-------------------------------------------

* ``BaseFab::atomicAdd`` and related helpers allow multiple tiles to accumulate
  into the same FAB in parallel (:file:`Src/Base/AMReX_BaseFab.H`).  On the GPU
  path the implementation goes through ``amrex::Gpu::Atomic::AddNoRet`` while
  the host path uses ``#pragma omp atomic`` in
  :file:`Src/Base/AMReX_GpuAtomic.H`.  The interleaving of those atomic updates
  depends on the thread and block schedule, so repeated runs can round the
  total differently.
* ``FabArray`` communication utilities rely on the same atomics when unpacking
  overlapping data.  For example ``FabArray<BaseFab>::ParallelCopy`` calls
  ``Gpu::Atomic::AddNoRet`` while combining receive buffers in
  :file:`Src/Base/AMReX_FBI.H`.  When multiple senders contribute to the same
  cell, the order of those additions is not deterministic.
* Embedded-boundary redistribution routines, such as
  ``EB_StateRedistribute`` and ``EB_StateRedistUtils``, update conservative
  quantities with GPU atomics (:file:`Src/EB/AMReX_EB_StateRedistribute.cpp`,
  :file:`Src/EB/AMReX_EB_StateRedistUtils.cpp`).  The choice of CUDA/HIP warp
  scheduling determines the summation order.

AMR flux registers
------------------

* ``YAFluxRegister`` employs ``HostDevice::Atomic::Add`` while accumulating
  fine fluxes back onto coarse faces (:file:`Src/Boundary/AMReX_YAFluxRegister_3D_K.H`).
  The atomic updates allow multiple fine patches to update the same coarse cell
  concurrently, so the final ordering depends on GPU thread/block scheduling or
  OpenMP interleaving. ``EBFluxRegister`` inherits this implementation and
  introduces additional atomic accumulations when redistributing embedded
  boundary fluxes (:file:`Src/EB/AMReX_EBFluxRegister_3D_C.H`).
* The legacy ``FluxRegister`` implementation coarsens fine fluxes inside a
  single thread before writing to the coarse register, avoiding atomic
  operations (:file:`Src/AmrCore/AMReX_FluxReg_3D_C.H`).  Likewise,
  ``EdgeFluxRegister`` (used for constrained-transport electric fields) loops
  over all fine contributions within the thread that owns the coarse edge, so
  its updates are deterministic (:file:`Src/Boundary/AMReX_EdgeFluxRegister.cpp`).
  The Fortran-facing ``FlashFluxRegister`` follows the same pattern and keeps
  its reductions local to the thread that owns the coarse face
  (:file:`Src/F_Interfaces/AmrCore/AMReX_FlashFluxRegister.H`).

Particle to mesh deposition
---------------------------

* Particle interpolation and charge deposition kernels (e.g.
  ``amrex_deposit_cic``) perform ``Gpu::Atomic::AddNoRet`` into nodal or cell
  centered fields (:file:`Src/Particle/AMReX_Particle_mod_K.H`).  The host path
  falls back to ``BaseFab::atomicAdd``
  (:file:`Src/Particle/AMReX_ParticleContainerI.H`,
  :file:`Src/Particle/AMReX_ParticleMesh.H`).  As particles are processed in
  parallel, the sequence of adds depends on launch configuration or OpenMP
  scheduling.
* Higher-order particle-mesh operators reuse atomics as they stencil deposit;
  ``ParticleInterpolators::ParticleToMesh`` invokes ``Gpu::Atomic::AddNoRet``
  for each weighted contribution (:file:`Src/Particle/AMReX_ParticleInterpolators.H`).
  Communication stages reuse the same buffers, so overlapping deposits inherit
  the non-deterministic order.

Parallel reductions and norms
-----------------------------

* ``ReduceOps`` / ``ReduceData`` implement parallel norms and sums for
  FAB-based data (:file:`Src/Base/AMReX_Reduce.H`,
  :file:`Src/Base/AMReX_ParReduce.H`).  Each reducer combines per-thread and
  per-block partial results in an order that depends on the execution backend
  (CUDA blocks, HIP wavefronts, or OpenMP threads).  Calls such as
  ``MultiFab::sum``, ``MultiFab::norm0`` and particle reductions inherit that
  behavior.
* ``HostDevice::Atomic::Add`` is used internally when reductions fall back to a
  host implementation, again relying on ``#pragma omp atomic`` in
  :file:`Src/Base/AMReX_GpuAtomic.H`.

MPI collective reductions
-------------------------

* ``ParallelDescriptor::ReduceReal*`` functions use ``MPI_Allreduce`` or
  ``MPI_Reduce`` (:file:`Src/Base/AMReX_ParallelDescriptor.H`).  MPI libraries
  are free to reorder partial sums while building their reduction trees, and
  the order can change with communicator topology, process pinning, or MPI
  version.  Cross-rank sums, norms, and max/min operations therefore vary at
  the level of machine epsilon even when each rank performs a deterministic
  local reduction.

Practices for reproducibility
-----------------------------

* Fix the execution configuration: run with a single MPI rank or a fixed MPI
  decomposition, disable OpenMP when possible, and use deterministic GPU launch
  parameters to minimize variability.
* Prefer deterministic algorithms (e.g., serial accumulation or Kahan/compensated
  summation) in application code where bitwise reproducibility is required.
* Capture reference answers with the same build options, compiler, and hardware
  whenever regression tests depend on floating-point tolerances.
* Document the expected tolerance in tests and diagnostics; accept that atomic
  updates can differ by several ulps unless a reproducible algorithm is used.
