.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


.. _sec:gpu:build:

Building GPU Support
====================

Building with GNU Make
----------------------

To build AMReX with GPU support, add ``USE_CUDA=TRUE`` to the 
``GNUmakefile`` or as a command line argument.  AMReX does not
require OpenACC or CUDA Fortran, but application codes can use them
if they are supported by the compiler.  For OpenACC support, add
``USE_ACC=TRUE``.  PGI, Cray and GNU compilers support OpenACC.  Thus, 
for OpenACC, you must use ``COMP=pgi``, ``COMP=cray`` or ``COMP=gnu``.
Only IBM and PGI support CUDA Fortran, which is also built when
``USE_CUDA=TRUE``.  OpenMP is currently not supported with CUDA.

Compiling AMReX with CUDA requires compiling the code through NVIDIA's 
CUDA compiler driver in addition to the standard compiler.  This driver
is called ``nvcc`` and it requires a host compiler to work through. 
The default host compiler for NVCC is GCC even if ``COMP`` is set to 
a different compiler.  One can change this by setting ``NVCC_HOST_COMP``.
For example, ``COMP=pgi`` alone will compile C/C++ codes with NVCC/GCC
and Fortran codes with PGI, and link with PGI.  Using ``COMP=pgi`` and 
``NVCC_HOST_COMP=pgi`` will compile C/C++ codes with PGI and NVCC/PGI.

You can use ``Tutorials/Basic/HelloWorld_C`` to test your programming
environment.  Building with:

.. highlight:: console

::

   make COMP=gnu USE_CUDA=TRUE

should produce an executable named ``main3d.gnu.DEBUG.CUDA.ex``.  You
can run it and that will generate results like:

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

Most GPU related classes and functions are in ``namespace Gpu``,
which is inside ``namespace amrex``. For example, the GPU configuration
class ``Device`` can be referenced to at ``amrex::Gpu::Device``. Other
important objects in the Gpu namespace include ``Gpu::AsyncFab``. 

For portability, AMReX defines some macros for CUDA function qualifiers
and they should be preferred to allow execution with ``USE_CUDA=FALSE``.
These include:

.. highlight:: c++

::

   #define AMREX_GPU_HOST        __host__
   #define AMREX_GPU_DEVICE      __device__
   #define AMREX_GPU_GLOBAL      __global__
   #define AMREX_GPU_HOST_DEVICE __host__ __device__

Note that when AMReX is not built with CUDA, these macros expand to
empty space.

When AMReX is compiled with ``USE_CUDA=TRUE``, the preprocessor 
macros ``AMREX_USE_CUDA`` and ``AMREX_USE_GPU`` are defined for 
conditional programming.  For PGI and IBM compilers, 
``AMREX_USE_CUDA_FORTRAN`` is also defined, as well as
``-DAMREX_CUDA_FORT_GLOBAL='attributes(global)'``,
``-DAMREX_CUDA_FORT_DEVICE='attributes(device)'``, and
``-DAMREX_CUDA_FORT_HOST='attributes(host)'`` so that CUDA Fortran
functions can be properly labelled.  When AMReX is compiled with
``USE_ACC=TRUE``, ``AMREX_USE_ACC`` is defined.

In addition to AMReX's preprocessor macros, CUDA provides the 
``__CUDA_ARCH__`` macro which is only defined when in device code.
``__CUDA_ARCH__`` should be used when a ``__host__ __device__``
function requires separate code for the CPU and GPU implementations.

.. ===================================================================

.. _sec:gpu:memory:

Memory Allocation
=================

To provide portability and improve memory allocation performance,
AMReX provides a number of memory pools.  When compiled without
CUDA, all :cpp:`Arena` s implement standard :cpp:`new` and 
:cpp:`delete` operators. Without CUDA, the :cpp:`Arena` s each 
allocate with a specific type of GPU memory:

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

The Arena object returned by these calls provides access
to two functions:

.. highlight:: c++

::

   void* alloc (std::size_t sz);
   void free (void* p);

:cpp:`The_Arena()` is used for memory allocation of data in
:cpp:`BaseFab`.  Therefore the data in a :cpp:`MultiFab` is placed in
unified memory and is accessible from both CPU host and GPU device.  
This allows application codes to develop their GPU capability
gradually.  :cpp:`The_Managed_Arena()` is a separate pool of
unified memory, that is distinguished from :cpp:`The_Arena()` for 
performance reasons.  If you want to print out the current memory usage
of the Arenas, you can call :cpp:`amrex::Arena::PrintUsage()`.

.. ===================================================================

.. _sec:gpu:classes:

GPU Safe Classes and Functions
==============================

AMReX GPU work takes place inside of MFIter and particle loops. 
Therefore, there are two ways classes and functions have been modified 
to interact with the GPU: 

1) A small number of functions used within these loops are labelled using
``AMREX_GPU_DEVICE`` and can be called on the device. This includes member 
functions, such as :cpp:`IntVect::type()`, as well as non-member functions,
such as :cpp:`amrex::min` and `amrex::max`. In specialized cases,
classes are labeled such that the object can be constructed, destructed 
and its functions can be implemented on the device, including ``IntVect``.

2) Functions that contain MFIter or particle loops have been rewritten
to contain device launches. For example, the :cpp:`FillBoundary`
function cannot be called from device code, but calling it from
CPU will launch GPU kernels if AMReX is compiled with GPU support. 

Most functions and objects have not been given a device version due to
CUDA restrictions and to ensure applications maintain an internally consistent
GPU strategy.

In this section, we discuss some examples of AMReX device classes and functions 
that are important for programming GPUs.

GpuArray
--------
:cpp:`std::array` is used throughout AMReX, however its functions are not defined
in device code. GpuArray is AMReX's built-in alternative. It uses a C array that
can be passed to the device by value and has device functions for the :cpp:`[]`
operator, :cpp:`size()` and :cpp:`data()` that returns a pointer to the C array.
GpuArray can be used whenever a fixed array of data needs to be passed to the GPU.

Additional functions have been created to return GpuArray instead of `std::array`, 
including `Box::CellSizeArray()`, `Box::InvCellSizeArray()` and `Box::length3d()`.

Various Vectors
---------------
HERE HERE HERE HERE

Get info from Particle section. (Mention primarily used for particles for now).
Variable size vectors. Resize ONLY on the CPU. Accessed to fixed size data by
passing pointer to data to device. Based on type, be careful when you do things
to the vector. (Host, Device & Managed / how to use.)

Box, IntVect and IndexType
--------------------------

In AMReX, :cpp:`Box`, :cpp:`IntVect` and :cpp:`IndexType` 
are classes for representing indices.  These classes and most of 
their member functions, including constructors and destructors,
have both host and device versions.  They can be used in
device code.  

Geometry
--------

AMReX's :cpp:`Geometry` class is not a GPU safe class.  However, we often need
to use geometric information such as cell size and physical coordinates
in GPU kernels.  To utilize :cpp:`Geometry` on the GPUs, the data is copied
into a GPU safe class that can be passed by value to GPU kernels. This class 
is called :cpp:`GeometryData`, which is created by calling 
:cpp:`Geometry::data()`. The accessor functions of :cpp:`GeometryData` are
identical to :cpp:`Geometry`.

One limitation of this strategy is that :cpp:`Geometry` cannot not be changed
on the device. :cpp:`GeometryData` holds a disposable copy of the data that 
does not synchronize with :cpp:`Geometry` after use. Therefore, only change 
:cpp:`Geometry` on the CPU, outside MFIter loops with GPU kernels to avoid
race conditions.

BaseFab, FArrayBox and IArrayBox
--------------------------------

:cpp:`BaseFab<T>`, :cpp:`IArrayBox` and :cpp:`FArrayBox`
have some GPU support.  They cannot be constructed in device code, but
a pointer to them can be passed to GPU kernels from CPU code.  Many
of their member functions can be used in device code as long as they
have been allocated in device memory. Some of the device member
functions include :cpp:`view`, :cpp:`dataPtr`, :cpp:`box`, 
:cpp:`nComp`, and :cpp:`setVal`.

All :cpp:`BaseFab<T>` objects in :cpp:`FabArray<FAB>` are allocated in
unifed memory, including :cpp:`MultiFab` and :cpp:`FArrayBox` which are
derived from :cpp:`BaseFab`. A :cpp:`BaseFab<T>` object created 
on the stack in CPU code cannot be used in GPU device code, because
the object is in CPU memory.  However, a :cpp:`BaseFab` created with
:cpp:`new` on the heap is GPU safe, because :cpp:`BaseFab` has its own
overloaded `:cpp:`operator new` that allocates memory from :cpp:`The_Arena()`,
a managed memory arena.  For example,

.. highlight:: c++

::

    // We are in CPU code

    FArrayBox cpu_fab(box,ncomp);
    // FArrayBox* p_cpu_fab = &(cpu_fab) cannot be used in GPU device code!

    FArrayBox* p_gpu_fab = new FArrayBox(box,ncomp);
    // FArrayBox* p_gpu_fab can be used in GPU device code.

amrex::min and amrex::max
-------------------------
ADD THIS IN!!!!

.. ===================================================================

.. _sec:gpu:launch:

Kernel Launch
=============

AMReX uses :cpp:`MFIter` to iterate over a :cpp:`MultiFab`.  Inside
the loop, we call functions to work on :cpp:`FArrayBox` objects (see
:ref:`sec:basics:mfiter`).  With GPU support enabled, kernels are launched
 inside the :cpp:`MFIter` loop.  A tutorial example can be found in
``Tutorials/GPU/Launch``.  How to launch using each possible GPU language
will be introduced:

Launching a C++ function
------------------------

The part launching a CUDA C++ kernel is shown below.

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

The code above works whether it is compiled for GPUs or CPUs. 
:cpp:`TilingIfNotGPU()` returns ``false`` in the GPU case to
turn off tiling so that GPU kernels have more compute work to do.
When tiling is off, :cpp:`tilebox()` returns the valid box of the 
:cpp:`FArrayBox` for that iteration.  The :cpp:`MultiFab::fabPtr` 
function takes :cpp:`MFIter` and returns a managed pointer that can be 
captured by an extended C++ lambda function the user defined in the 
``AMREX_LAUNCH_DEVICE_LAMBDA`` launch macro.

The launch macro takes three arguments.  The first argument is a 
:cpp:`Box` specifying the whole region of the kernel.  The second 
argument is the desired name of a :cpp:`Box` variable that
specifies the subregion each thread works on.  This subregion 
:cpp:`Box`  is defined with the specificed name inside the macro.  
The third  argument denotes the code block of the lambda function 
that will be ran in this launch. In this example, a GPU device 
function, ``plusone_cudacpp``, is called and is passed the captured 
variable ``fab``.  In CUDA, an extended lambda function can only 
capture by value, not reference.  That's why a pointer must be created
to the :cpp:`FArrayBox`.

Launching a FORTRAN function
----------------------------

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

Offloading work using OpenACC pragmas
-------------------------------------

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

Note that here:cpp:`MultiFab::operator[]` is used to get a reference
to :cpp:`FArrayBox` rather than :cpp:`MultiFab::fabPtr` to get
a pointer, as was done in the CUDA examples.  Using the reference is 
optimal for performance when passing pointers to kernels is not 
required, which includes CPU code sections and OpenACC implementations. 

Function ``plusone_acc`` is a CPU host function.  The reference from 
:cpp:`operator[]` is a reference to a :cpp:`FArrayBox` in host
memory even though the data pointer inside the object points to
unified memory.  This managed data pointer is retrieved with 
:cpp:`fabPtr`.  ``BL_TO_FORTRAN_ANYD`` expands to the individual
components of the :cpp:`FArrayBox`, including the :cpp:`Box` defining
its indicies, the number of components and the data pointer itself.
By passing the :cpp:`FArrayBox` via its required components, 
unnecessary data movement is minimized. 

The corresponding OpenACC labelled loop in ``plusone_acc`` is: 

.. highlight:: fortran 

::
    !pointer to fab managed data defined as "dat"

    !$acc kernels deviceptr(dat)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat(i,j,k) = dat(i,j,k) + 1.0_amrex_real
          end do
       end do
    end do
    !$acc end kernels


Since the data pointer passed to ``plusone_acc`` points to 
unified memory, OpenACC is told the data is available on the device
by using the ``deviceptr`` construct.

Launching over a particle loop
------------------------------

SIMILAR PARTICLE LAUNCH EXAMPLE HERE???
In the three examples above, we have shown how to offload to the GPU.
In the CUDA C++ and Fortran cases, the kernel is launched with a C++
template global function in AMReX (hidden in the launch macro),
whereas for the OpenACC example, it is done with pragma in Fortran.
More details on the examples can be found in the source codes at
``Tutorials/GPU/Launch/``.  We also refer the readers to Chapter
:ref:`Chap:Basics` for information about basic AMReX classes.

Kernel launch details
---------------------

CUDA kernel calls are asynchronous and they return before the kernel 
is finished on the GPU. So :cpp:`MFIter` finishes its iterations on
the CPU before the GPU finishes its work. To guarantee data coherence,
there is an implicit CUDA device synchronization (a CUDA barrier) in 
the destructor of :cpp:`MFIter`. This makes each iteration of 
:cpp:`MFIter` independent, but potentially concurrent.

CUDA supports multiple streams and kernels. The kernels in the 
same stream are executed sequentially, but kernels from different streams
may be run in parallel.  For each iteration of :cpp:`MFIter`, AMReX 
uses a different CUDA stream (up to 16 streams in total) for the 
kernels in that iteration. 

These  
allows each :cpp:`FArrayBox` inside of a :cpp:`MultiFab` to be implemented.

SECTION ON EXTENDED LAMBDA RESTRICTIONS

Launching kernels with the ``AMREX_LAUNCH_DEVICE_LAMBDA`` uses the CUDA
extended lamdba feature.  There are, however, some restrictions on
extended lambdas the user must be aware of.  For example, the enclosing
function for the extended lamdba must not have private or protected
access within its parent class, otherwise the code will not compile.
In that case, we have to change the function to a public member of its
parent class.  There is also a pitfall we *must* be aware.  If the
extended lambda accesses a member of the enclosing class, the lambda
function actually captures :cpp:`this` pointer by value and accesses
variable via :cpp:`this->`.  If the object is not accessible on GPU,
the code will not work as intended.  For example,

.. highlight:: c++

::

    class MyClass {
    public:
        Box bx;
        int m;
        void f () {
            AMREX_LAUNCH_DEVICE_LAMBDA (bx, tbx,
            {
                printf("m = %d\n", m);
            });
        }
    };

The function ``f`` in the code above will not work unless :cpp:`MyClass`
object is in unified memory.  To fix it, we can change it to

.. highlight:: c++

::

    class MyClass {
    public:
        Box bx;
        int m;
        void f () {
            int local_m = m;
            AMREX_LAUNCH_DEVICE_LAMBDA (bx, tbx,
            {
                printf("m = %d\n", local_m);
            });
        }
    };


.. _sec:gpu:example:

An Example of Migrating to GPU
==============================

The nature of GPU programming poses difficulties for a common patterns
like the one below.

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

There are several issues in migrating the code above to GPU that need to
be addressed.  Because functions ``f1`` and ``f2`` have different
work regions and there are data dependencies between the two, it is
difficult to put them into a single GPU kernel.  So, we will launch two
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
will compete for access to ``q``.  Moving the line into the body of
:cpp:`MFIter` loop will make ``q`` a variable local to each iteration,
but it has a new issue.  When do we delete :cpp:`q`?  To the CPU, the
resource of :cpp:`q` should be freed at the end of the scope, otherwise
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

AMReX provides some functions for performing reduction operations with on the GPU (e.g.,
:cpp:`MultiFab::sum`, :cpp:`MultiFab::max`, etc.). Function
templates :cpp:`amrex::ReduceSum`, :cpp:`amrex::ReduceMin` and
:cpp:`amrex::ReduceMax` can be used to implement your own reduction
functions for :cpp:`FabArray`\ s.

.. ===================================================================

Particle Support
================

.. _sec:gpu:particle:

AMReX's GPU particle support relies on Thrust, a parallel algorithms library maintained by
Nvidia. Thrust provides a GPU-capable vector container that is otherwise similar to the one
in the C++ Standard Template Library, along with associated sorting, searching, and prefix
summing operations. Along with Cuda's unified memory, thrust forms the basis of our GPU support
for particles. 

When compiled with ``USE_CUDA=TRUE``, AMReX places all its particle data in instances of
``thrust::device_vector`` that have been configured using a custom memory allocator that
wraps around ``cudaMallocManaged``. This means that the :cpp:`dataPtr` associated with particle
data can be passed in to GPU kernels, similar to way it would be passed in to a Fortran
subroutine in typical AMReX CPU code. As with the mesh data, these kernels can be launched
with a variety of approaches, including Cuda C / Fortran and OpenACC.

For portability, we have provided a set of Vector classes that wrap around the Thrust and
STL vectors. When ``USE_CUDA = FALSE``, these classes reduce to the normal :cpp:`amrex::Vector`.
When ``USE_CUDA = TRUE``, they have different meanings. :cpp:`Gpu::HostVector` is a wrapper
around :cpp:`thrust::host_vector`. :cpp:`Gpu::DeviceVector` is a wrapper around :cpp:`thrust::device_vector`,
while :cpp:`Gpu::ManagedDeviceVector` is a :cpp:`thrust::device_vector` that lives in managed memory. These classes are useful when you have certain stages of an algorithm that you know will always
execute on either the host or the device.

AMReX's :cpp:`Redistribute()`, which moves particles back to the proper grids after their positions
have changed, has been ported to work on the GPU as well. You can't call it from device code,
but you can call it on particles that reside on the device and it won't trigger any unified
memory traffic. As with :cpp:`MultiFab` data, the MPI portion of the particle redistribute is set
up to take advantange of the Cuda-aware MPI implementations available on platforms such as
ORNL's Summit and Summit-dev.

Performance Tips
================

.. _sec:gpu:performance:

It may be helpful to keep in mind these performance tips when working with AMReX for GPUs:

*
*
*

.. ===================================================================

Limitations
===========

.. _sec:gpu:limits:

GPU support in AMReX is still under development.  There are some known
limitations:

- By default, AMReX assumes the MPI library used is GPU aware.  The
  communication buffers given to MPI functions are allocated in device
  memory.

- OpenMP is currently not compatible with building AMReX with CUDA. 
  ``USE_CUDA=TRUE`` and ``USE_OMP=TRUE`` will fail to compile.

- CMake is not yet supported for building AMReX GPU support.

- Many multi-level functions in AMReX have not been ported to GPUs.

- Linear solvers have not been ported to GPUs.

- Embedded boundary capability has not been ported to GPUs.

- The Fortran interface of AMReX does not currently have GPU support.
