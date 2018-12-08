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

See the source codes in ``Tutorials/GPU/Launch/`` for more details on
the kernels.

.. streams and mfiter

.. inLaunchRegion

.. macro, cuda fortran and openacc

.. refer back to Basic for calling C++ functions on FArrayBox

.. ===================================================================

.. _sec:gpu:safeclasses:

GPU Safe Classes
================

.. ===================================================================

.. _sec:gpu:assertion:

Assertion and Error Check
=========================

.. AMREX_GPU_SAFE_CALL(), AMREX_GPU_ERROR_CHECK();

.. ===================================================================

.. _sec:gpu:reduction:

Reduction
=========

.. ===================================================================

.. _sec:gpu:particle:

Particle
========

.. ===================================================================

.. _sec::gpu:mpi:

CUDA Aware MPI
==============

.. ===================================================================

.. _sec:gpu:limits:

Pitfalls and Limitations
========================

.. At most one gpu per mpi rank.  OpenMP, AMR development are underway
.. cmake and build as library
