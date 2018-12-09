.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


.. _sec:gpu:build:

Building GPU Support
====================

Building with GNU Make
----------------------

To build AMReX with GPU support, you add ``USE_CUDA=TRUE`` to your
make file or as an command line argument.  AMReX itself does not
require OpenACC and CUDA Fortran, but application codes can use them
if they are supported by the compiler.  For OpenACC support, you add
``USE_ACC=TRUE``.  Only IBM and PGI support CUDA Fortran.  Thus for
CUDA Fortran, you must use ``COMP=pgi`` or ``COMP=ibm``.  OpenMP is
currently not supported when ``USE_CUDA=TRUE``.  The default host
compiler for NVCC is GCC even if ``COMP`` is set to a different
compiler.  One can change this by setting ``NVCC_HOST_COMP`` to say
``pgi``.  For example, ``COMP=pgi`` alone will compile C/C++ codes
with NVCC/GCC and Fortran codes with PGI, and link with PGI.  Using
``COMP=pgi`` and ``NVCC_HOST_COMP=pgi`` will compile C/C++ codes with
NVCC/PGI.

You can use ``Tutorials/Basic/HelloWorld_C`` to test your programming
environment.  Typing

.. highlight:: console

::

   make COMP=gnu USE_CUDA=TRUE

should produce an executable named ``main3d.gnu.DEBUG.CUDA.ex``.  You
can run it and that will generate results like below

.. highlight:: console

::

   $ ./main3d.gnu.DEBUG.CUDA.ex 
   CUDA initialized with 1 GPU
   AMReX (18.12-95-gf265b537f479-dirty) initialized
   Hello world from AMReX version 18.12-95-gf265b537f479-dirty
   [The         Arena] space (kilobyte): 8192
   [The  Device Arena] space (kilobyte): 8192
   [The Managed Arena] space (kilobyte): 8192
   [The  Pinned Arena] space (kilobyte): 8192
   AMReX (18.12-95-gf265b537f479-dirty) finalized

.. ===================================================================

.. _sec:gpu:namespace:

Gpu Namespace and Macros
========================

GPU related classes and functions are usually in ``namespace Gpu``,
which is inside ``namespace amrex``.  For portability, AMReX defines
some macros for CUDA function qualifiers and they should be preferred
to hardwired ``__*__``.  These include

.. highlight:: c++

::

   #define AMREX_GPU_HOST        __host__
   #define AMREX_GPU_DEVICE      __device__
   #define AMREX_GPU_GLOBAL      __global__
   #define AMREX_GPU_HOST_DEVICE __host__ __device__

Note that when AMReX is not built with CUDA, these macros expand to
empty space.

When AMReX is compiled with ``USE_CUDA=TRUE``, we pass
``-DAMREX_USE_CUDA`` and ``-DAMREX_USE_GPU`` to the compiler so that
these macros can be used for preprocessing.  For PGI and IBM
compilers, we also pass ``-DAMREX_USE_CUDA_FORTRAN``,
``-DAMREX_CUDA_FORT_GLOBAL='attributes(global)'``,
``-DAMREX_CUDA_FORT_DEVICE='attributes(device)'``, and
``-DAMREX_CUDA_FORT_HOST='attributes(host)'`` so that CUDA Fortran
functions can be properly declared.  When AMReX is compiled with
``USE_ACC=TRUE``, we pass ``-DAMREX_USE_ACC`` to the compiler.

.. ===================================================================

.. _sec:gpu:memory:

Memory Allocation
=================

To provide portability and improve memory allocation performance,
AMReX provides a number of memory pools.

.. raw:: latex

    \begin{center}

.. _tab:gpu:arena:

.. table:: Memory Arenas

    +---------------------+------------------+
    | Arena               |    Memory Type   |
    +=====================+==================+
    | The_Arena()         |  unified memory  | 
    +---------------------+------------------+
    | The_Device_Arena()  |  device memory   | 
    +---------------------+------------------+
    | The_Managed_Arena() |  unified memory  | 
    +---------------------+------------------+
    | The_Pinned_Arena()  |  pinned memory   | 
    +---------------------+------------------+

.. raw:: latex

    \end{center}

:cpp:`Arena` object returned by these arena functions provides

.. highlight:: c++

::

   void* alloc (std::size_t sz);
   void free (void* p);

:cpp:`The_Arena()` is used for memory allocation of data in
:cpp:`BaseFab`.  Therefore the data in a :cpp:`MultiFab` are in
unified memory and they are accessible from both CPU host and GPU
device.  This allows application codes to develop their GPU capability
gradually.  :cpp:`The_Managed_Arena()` is also a memory pool of
unified memory, but it is separated from :cpp:`The_Arena()` for
performance reason.  If you want to print out the current memory usage
by the arenas, you can call :cpp:`amrex::Arena::PrintUsage()`.

.. ===================================================================

.. _sec:gpu:classes:

GPU Safe Classes and Functions
==============================

Some AMReX classes and functions can be used in device functions, but
most do not.  Note that this is about using them on GPU, not about
whether they use GPU.  For example, :cpp:`MultiFab::Copy` function
cannot be called from GPU, but calling it from CPU will use GPU if
AMReX is compiled with GPU support.  In this section, we discuss a
few classes and functions that are useful for programming GPU.

In the basic AMReX classes, :cpp:`Box`, :cpp:`IntVect` and
:cpp:`IndexType` are classes for representing indices.  These classes
and most of their member functions, including constructors and
destructors, have both host and device versions.  They can be used in
device code.  :cpp:`BaseFab<T>`, :cpp:`IArrayBox` and :cpp:`FArrayBox`
have some GPU support.  They cannot be constructed in device code, but
a pointer to them can be passed to GPU kernels from CPU code.  Many
member functions of them (e.g., :cpp:`view`, :cpp:`dataPtr`,
:cpp:`box`, :cpp:`nComp`, and :cpp:`setVal`) can be used in device
code, if the pointer points to unified memory.  All :cpp:`BaseFab<T>`
objects (including :cpp:`FArrayBox` derived from :cpp:`BaseFab`) in
:cpp:`FabArray<FAB>` (including :cpp:`MultiFab`) are allocated in
unified memory.  A :cpp:`BaseFab<T>` object created on the stack in
CPU code cannot be used in GPU device code, because the object is in
CPU memory.  However, a :cpp:`BaseFab` created with :cpp:`new` on the
heap is GPU safe, because :cpp:`BaseFab` has its own overloaded
`:cpp:`operator new` that allocates memory from :cpp:`The_Arena()`, a
managed memory arena.  For example,

.. highlight:: c++

::

    // We are in CPU code

    FArrayBox cpu_fab(box,ncomp);
    // FArrayBox* p_cpu_fab = &(cpu_fab) cannot be used in GPU device code!

   FArrayBox* p_gpu_fab = new FArrayBox(box,ncomp);
   // FArrayBox* p_gpu_fab can be used in GPU device code.

:cpp:`Geometry` class is not a GPU safe class.  However, we often need
to use geometry information such as cell size and physical coordinates
in GPU kernels.  What we can do is extract its data into a GPU safe
class :cpp:`GeometryData` with :cpp:`Geometry::data` function and pass
it by value to GPU kernels.

.. ===================================================================

.. _sec:gpu:launch:

Kernel Launch
=============

AMReX uses :cpp:`MFIter` to iterate over a :cpp:`MultiFab`.  Inside
the loop, we call functions to work on :cpp:`FArrayBox` objects (see
:ref:`sec:basics:mfiter`).  With GPU, we launch kernels inside
:cpp:`MFIter` loop.  A tutorial example can be found in
``Tutorials/GPU/Launch``.  The part launching a CUDA C++ kernel is
shown below.

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox* fab = mf.fabPtr(mfi);
        AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        {
            plusone_cudacpp(tbx, *fab);
        });
    }

The code above works whether it is compiled for GPU or CPU.  We use
:cpp:`TilingIfNotGPU()` that returns ``false`` in the case of GPU to
turn off tiling so that GPU kernels have more works.  When tiling is
off, :cpp:`tilebox()` returns the valid box of the :cpp:`FArrayBox`
for that iteration.  :cpp:`MultiFab::fabPtr` function takes
:cpp:`MFIter` and returns a managed pointer that is subsequently
captured by an extended C++ lambda function produced by the
``AMREX_LAUNCH_DEVICE_LAMBDA`` macro.  The launch macro usually takes
three arguments.  In this example, the first argument is a :cpp:`Box`
specifying the whole region of the kernel.  The second argument is a
:cpp:`Box` variable that specifies the region a thread works on.  Note
that the second argument is the name of a local variable to the thread
defined by the macro, not an existing variable.  The third argument is
a block of codes delimited by a pair of curly braces.  In that block,
we call a GPU device function ``plusone_cudacpp`` with captured
variable ``fab``.  In CUDA, an extended lambda function can only
capture by value, not reference.  That's why we capture a pointer to
:cpp:`FArrayBox`.

We can also call CUDA Fortran device functions in the code block for
the launch macro like below.

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox* fab = mf.fabPtr(mfi);
        AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        {
            plusone_cudafort(BL_TO_FORTRAN_BOX(tbx),
                             BL_TO_FORTRAN_ANYD(*fab));
        });
    }

Because :cpp:`Box` and :cpp:`FArrayBox` are C++ classes not understood by
Fortran, we use some helper macros to pass them as Fortran data types
(see :ref:`sec:basics:fortran`).

The tutorial at ``Tutorials/GPU/Launch`` also shows an example of
using OpenACC in Fortran.  We call a Fortran function and in that
function we use OpenACC to offload work to GPU.

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        plusone_acc(BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_ANYD(fab));
    }

Note that here we use :cpp:`MultiFab::operator[]` to get a reference
to :cpp:`FArrayBox` as what we usually do for CPU codes, rather using
:cpp:`MultiFab::fabPtr` to get a pointer for the CUDA examples we just
showed above.  The reason for this is performance.  Function
``plusone_acc`` is a CPU host function.  The reference we get from
:cpp:`operator[]` is a reference to a :cpp:`FArrayBox` in host memory
even though its data pointer inside the object points to unified
memory, whereas the pointer we get from :cpp:`fatPtr` is a manged
memory pointer.  Because ``BL_TO_FORTRAN_ANYD`` in this case expands
to the CPU version of some :cpp:`FArrayBox` member functions (unlike
GPU functions in the CUDA Fortran example above), having the metadata
(i.e., :cpp:`Box`, the number of components and the data pointer
itself) can minimize unnecessary data movement.  Since the data
pointer passed to ``plusone_acc`` as Fortran array by the
``BL_TO_FORTRAN_AND`` macro points to unified memory, we can take
advantage of that by declaring it as OpenACC ``deviceptr``.

In the three examples above, we have shown how to launch GPU kernels.
In the CUDA C++ and Fortran cases, the kernel is launched with a C++
template global function of AMReX (hidden in the launch macro),
whereas for the OpenACC example, it is done with pragma in Fortran.
More details on the examples can be found in the source codes at
``Tutorials/GPU/Launch/``.  We also refer the readers to Chapter
:ref:`Chap:Basics` for information about basic AMReX classes.

CUDA supports multiple streams and kernels.  For each iteration of
:cpp:`MFIter`, AMReX uses a different stream (up to 16 streams in
total) for the kernels in that iteration.  The kernels in the same
stream are executed sequentially, but kernels from different streams
may be run in parallel.  Note that, to CPU, CUDA kernel calls are
asynchronous and they return before the kernel is finished on GPU.  So
:cpp:`MFIter` finishes its iterations on CPU before GPU finishes its
work.  To guarantee data coherence, there is an implicit CUDA device
synchronization in the destructor of :cpp:`MFIter`.

.. _sec:gpu:example:

An Example of Migrating to GPU
==============================

The nature of GPU programming poses difficulty for a common pattern
like below.

.. highlight:: c++

::

   // Given MultiFab uin and uout
   FArrayBox q;
   for (MFIter mfi(uin); mfi.isValid(); ++mfi)
   {
       const Box& vbx = mfi.validbox();
       const Box& gbx = amrex::grow(vbx,1);
       q.resize(gbx);

       // Do some work with uin[mfi] as input and q as output.
       // The output region is gbx;
       f1(gbx, q, uin[mfi]);

       // Then do more work with q as input and uout[mfi] as output.
       // The output region is vbx.
       f2(vbx, uout[mfi], q);
   }

There are several issues of migrating the code above to GPU needed to
been addressed.  Because functions ``f1`` and ``f2`` have different
work regions and there are data dependency between the two, it is
difficult to put them into a single GPU kernel.  So we will launch two
separate kernels.

As we have discussed in Section :ref:`sec:gpu:classes`, all
:cpp:`FArrayBox`\ es in the two :cpp:`MultiFab`\ s are in unified
memory.  But :cpp:`FArrayBox q` is in host memory.  Changing it to

.. highlight:: c++

::

    FArrayBox* q = new FArrayBox;

does not solve the problem completely because GPU kernel calls are
asynchronous from CPU's point of view.  Therefore there is a race
condition that GPU kernels in different iterations of :cpp:`MFIter`
will compete the access to ``q``.  Moving the line into the body of
:cpp:`MFIter` loop will make ``q`` a variable local to each iteration,
but it has a new issue.  When do we delete :cpp:`q`?  To CPU, the
resource of :cpp:`q` should be freed at the end of the scope otherwise
there will be a memory leak.  But at the end of the CPU scope, GPU
kernels might still need it.

One way to fix this is put the temporary :cpp:`FArrayBox` objects in a
:cpp:`MultiFab`.  Another way is to use :cpp:`Gpu:AsyncFab` designed
for this kind of situation.  In the code below, we show how it is used
and how kernels are launched.

.. highlight:: c++

::

   for (MFIter mfi(uin); mfi.isValid(); ++mfi)
   {
       const Box& vbx = mfi.validbox();
       const Box& gbx = amrex::grow(vbx,1);
       Gpu::AsyncFab q(gbx);
       FArrayBox const* uinfab  = uin.fabPtr();
       FArrayBox      * uoutfab = uout.fabPtr();

       AMREX_LAUNCH_DEVICE_LAMBDA ( gbx, tbx,
       {
           f1(tbx, q.fab(), *uinfab);
       };

       AMrEX_LAMBDA_DEVICE_LAMBDA ( vbx, tbx,
       {
           f2(tbx, *uoutfab, q.fab());
       });
   }

.. ===================================================================

.. _sec:gpu:assertion:

Assertion and Error Check
=========================

To help debugging, we often use :cpp:`amrex::Assert` and
:cpp:`amrex::Abort`.  These functions are GPU safe and can be used in
GPU kernels.  In CPU code, we can also call
:cpp:`AMREX_GPU_ERROR_CHECK()` to check if there are any GPU errors at
that point.  Because of asynchronicity, even if GPU kernels launched
before that point contain a bug that will result in a CUDA error, the
error may not be encountered when :cpp:`AMREX_GPU_ERROR_CHECK()` is
called.  We can use :cpp:`Gpu::Device::synchronize()` or
:cpp:`Gpu::Device::streamSynchroniz()` to synchronize the device or
the CUDA stream, respectively.

.. ===================================================================

.. _sec:gpu:reduction:

Reduction
=========

AMReX provides some functions for performing reduction with GPU (e.g.,
:cpp:`MultiFab::sum`, :cpp:`MultiFab::`max`, etc.).   Function
templates :cpp:`amrex::ReduceSum`, :cpp:`amrex::ReduceMin` and
:cpp:`amrex::ReduceMax` can be used to implement your own reduction
functions for :cpp:`FabArray`\ s.

.. ===================================================================

.. .. _sec:gpu:particle:

.. Particle
.. ========

.. ===================================================================

.. _sec:gpu:limits:

Limitations
===========

GPU support in AMReX is still under development.  There are some know
limitations.

- By default, AMReX assumes the MPI library used is GPU aware.  The
  communication buffers given to MPI functions are allocated in device
  memory.

- OpenMP is currently not compatible with building AMReX
  with ``USE_CUDA=TRUE``.

- CMake is not yet supported for building AMReX GPU support.

- Many multi-level functions in AMReX have not been ported to GPU.

- Linear solvers have not been ported to GPU.

- Embedded boundary capability has not been ported to GPU.

- Fortran interface of AMReX does not currently have GPU support.
