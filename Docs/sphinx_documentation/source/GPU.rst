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

Building with CMake
-------------------

CMake is currently unavailable when building with GPUs.

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

#. A small number of functions used within these loops are labelled using
``AMREX_GPU_DEVICE`` and can be called on the device. This includes member 
functions, such as :cpp:`IntVect::type()`, as well as non-member functions,
such as :cpp:`amrex::min` and `amrex::max`. In specialized cases,
classes are labeled such that the object can be constructed, destructed 
and its functions can be implemented on the device, including ``IntVect``.

#. Functions that contain MFIter or particle loops have been rewritten
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

ManagedVector
-------------

AMReX also provides a dynamic memory allocation object for GPU managed memory:
:cpp:`Gpu::ManagedVector`.  This class behaves identically to an
:cpp:`amrex::Vector`, (see :ref:`sec:basics:vecandarr`), except the vector's 
allocator has been changed to allocate and deallocate its data in CUDA
managed memory whenever ``USE_CUDA=TRUE``.

While the data is managed and available on GPUs, the member functions of
:cpp:`Gpu::ManagedVector` are not. To use the data on the GPU, it is
necessary to pass the underlying data pointer to the GPU. The managed data
pointer can be accessed using the :cpp:`data()` member function. 

Be aware: resizing of dynamically allocated memory on the GPU is unsupported.
All resizing of the vector should be done on the CPU, in a manner that avoids
race conditions with concurrent GPU kernels. 

CUDA's Thrust Vectors 
---------------------

CUDA's Thrust library can also be used to manage dynamically sized data sets.
However, if Thrust is used directly in AMReX code, it will be unable to compile
for cases when ``USE_CUDA=FALSE``.  To alleviate this issue, 
:cpp:`thrust::host_vector` and :cpp:`thrust::device_vector` have been wrapped
into the AMReX classes :cpp:`Gpu::HostVector` and :cpp:`Gpu::DeviceVector`,
When ``USE_CUDA=FALSE``, these classes revert to AMReX's Vector class. When 
``USE_CUDA=TRUE``, these classes become the corresponding Thrust vector.

Just like with Thrust vectors, :cpp:`HostVector` and :cpp:`DeviceVector` cannot 
be directly used on the device. For convenience, the :cpp:`dataPtr()` member
function has been altered to implement :cpp:`thrust::raw_pointer_cast` and 
return the raw data pointer which can be used to access the vector's underlying
data on the GPU.

It has proven useful to have a version of Thrust's :cpp:`device_vector` 
that uses CUDA managed memory. This is provided by :cpp:`Gpu::ManagedDeviceVector`. 

:cpp:`thrust::copy` is also commonly used in AMReX applications. It can be
implemented portably using :cpp:`Gpu::thrust_copy`. 

amrex::min and amrex::max
-------------------------

GPU versions of ``std::min`` and ``std::max`` are not provided in CUDA.
So, AMReX provides a templated :cpp:`min` and :cpp:`max` with host and 
device versions to allow functionality on GPUs. Invoke the explicitly 
namespaced :cpp:`amrex::min(A, B)` or :cpp:`amrex::max(x, y)` to use the 
GPU safe implementations. 


MultiFab Reductions
-------------------

AMReX provides functions for performing standard reduction operations on 
:cpp:`Multifabs`, including :cpp:`MultiFab::sum` and :cpp:`MultiFab::max`.
When ``USE_CUDA=TRUE``, these functions automatically implement the 
corresponding reductions on GPUs in a highly efficient manner.

Function templates :cpp:`amrex::ReduceSum`, :cpp:`amrex::ReduceMin` and
:cpp:`amrex::ReduceMax` can be used to implement user-defined reduction
functions over :cpp:`MultiFab`\ s. These same templates are implemented 
in the :cpp:`Multifab` functions, so they can be used as a reference to
build a custom reduction. For example, the :cpp:`Multifab:dot` 
implementation is reproduced here:

.. highlight:: c++

::

    Real sm = amrex::ReduceSum(x, y, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
    {
        return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp);
    });

    if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;

:cpp:`amrex::ReduceSum` takes two :cpp:`MultiFabs`, ``x`` and ``y`` and
returns the sum of the value returned from the given lambda function.
In this case, :cpp:`BaseFab::dot` is returned, yielding a sum of the
dot product of each local pair of :cpp:`BaseFab` s. Finally, 
:cpp:`ParallelAllReduce` is used to sum the dot products across all
MPI ranks and return the total dot product of the two :cpp:`MultiFabs`.

To implement a different reduction, replace the code block inside the
lambda function with the operation that should be applied, being sure
to return the value to be summed, minimized, or maximized.  The reduction
templates have a few different interfaces to accomodate a variety of 
reductions.  The :cpp:`amrex::ReduceSum` reduction template has varieties
that take either one, two or three ::cpp:`MultiFab` s. 
:cpp:`amrex::ReduceMin` and :cpp:`amrex::ReduceMax` can take either one
or two.

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
:cpp:`Geometry` on the CPU and outside of MFIter loops with GPU kernels to
avoid race conditions.

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
overloaded `:cpp:`operator new` that allocates memory from
:cpp:`The_Arena()`, a managed memory arena.  For example,

.. highlight:: c++

::

    // We are in CPU code

    FArrayBox cpu_fab(box,ncomp);
    // FArrayBox* p_cpu_fab = &(cpu_fab) cannot be used in GPU device code!

    FArrayBox* p_gpu_fab = new FArrayBox(box,ncomp);
    // FArrayBox* p_gpu_fab can be used in GPU device code.

.. ===================================================================

.. _sec:gpu:launch:

Kernel Launch
=============

In this section, how to offload work to the GPU will be demonstrated.
In the CUDA C++ and Fortran cases, the kernel is launched with a C++
template global function in AMReX (hidden in an AMReX launch macro),
whereas for the OpenACC example, it is done with pragmas in Fortran.
More details on the examples can be found in the source codes at
``Tutorials/GPU/Launch/``.  We also refer the readers to Chapter
:ref:`Chap:Basics` for information about basic AMReX classes.

It is also recommended to write your primary floating point operation 
kernels in C++ using AMReX's :cpp:`FabView` object syntax. It converts 
the 1D array into a simple to understand 3D loop structure, similar
to Fortran while maintaining performance. The details can be found in
:ref:`C++ Kernel <sec:basics:cppkernel>`.

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

    !dat = pointer to fab's managed data 

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

Launching an generic loop
-------------------------

The previous examples showed how the ``AMREX_LAUNCH_DEVICE_LAMBDA``
macro can be used to launch threads that work across cells in a
:cpp:`FArrayBox` inside of an :cpp:`MFIter` loop. However, the macro
is also capable of launching for a specified number of iterations
that will be split across GPU threads. For example, launching over
the number of elements in a vector is an example given in 
``Tutorials/GPU/Launch`` and is reproduced here: 

.. highlight:: c++

::

    {
        int size = 100;
        amrex::Gpu::ManagedVector<int> ones(size, 0);
        const auto data = ones.dataPtr();
        AMREX_LAUNCH_DEVICE_LAMBDA(size, iter,
        {
            data[iter] = data[iter] + 1;
        });

        Gpu::Device::synchronize();
    }


A :cpp:`ManagedVector` is created that can be worked on from 
both the CPU and GPU.  A copy of the underlying data pointer is
created, so it can be captured and passed into the lambda function.
The ``AMREX_LAUNCH_DEVICE_LAMBDA`` launch macro is used and again 
it takes three arguments: the size of the vector as a :cpp:`long`,
a name for the loop iterator and a lambda function to perform on
each iteration.

This form of the macro can be used to create any standard, singly 
incremented loop on a GPU, not just over a vector.  This macro
can also work on any contiguous subsets of a data set by passing
a pointer that points to the beginning of the subset and the size
of the subset.

Users must also be aware that there is no automatic device 
synchronize, as these launches do not have to occur inside of an 
:cpp:`MFIter` loop.  Add :cpp:`Gpu::Device::synchronize()` where
necessary to ensure GPU and CPU resources to do alter the same
data simultaneously. 

Kernel launch details
---------------------

CUDA kernel calls are asynchronous and they return before the kernel 
is finished on the GPU. So :cpp:`MFIter` finishes its iterations on
the CPU before the GPU finishes its work.  To guarantee consistency,
there is an implicit CUDA device synchronization (a CUDA barrier) in 
the destructor of :cpp:`MFIter`.  This ensures that all GPU work
inside of an :cpp:`MFIter` loop will complete before code outside of
the loop is executed.

CUDA supports multiple streams and kernels. Kernels launched in the 
same stream are executed sequentially, but different streams of kernel
launches may be run in parallel.  For each iteration of :cpp:`MFIter`, 
AMReX uses a different CUDA stream (up to 16 streams in total) for the 
kernels in that iteration.  This allows each iteration of an :cpp:`MFIter`
loop to run indepenently and maximize the use of GPU resources while
writting clean, readable :cpp:`MFIter` loops.

Launching kernels with the ``AMREX_LAUNCH_DEVICE_LAMBDA`` uses the CUDA
extended lamdba feature.  Extended lambdas have some restrictions the user
must understand.  For example, the function enclosing the extended
lamdba must not have private or protected access within its parent class, 
otherwise the code will not compile.  This can be fixed by changing the
access of the enclosing function to public. 

Another pitfall that *must* be considered: if the
extended lambda accesses a member of the enclosing class, the lambda
function actually captures :cpp:`this` pointer by value and accesses
variable via :cpp:`this->`.  If the object is not accessible on GPU,
the code will not work as intended.  For example,

.. highlight:: c++

::

    class MyClass {
    public:
        Box bx;
        int m;                           // Unmanaged integer created on the host.
        void f () { 
            AMREX_LAUNCH_DEVICE_LAMBDA (bx, tbx,
            {
                printf("m = %d\n", m);   // Failed attempt to use m on the GPU.
            });
        }
    };

The function ``f`` in the code above will not work unless :cpp:`MyClass`
object is in unified memory.  If it is undesirable to put the class into
unified memory, a local copy of the information can be created for the
lambda to capture. For example:

.. highlight:: c++

::

    class MyClass {
    public:
        Box bx;
        int m;
        void f () {
            int local_m = m;                  // Local temporary copy of m.
            AMREX_LAUNCH_DEVICE_LAMBDA (bx, tbx,
            {
                printf("m = %d\n", local_m);  // Lambda captures local_m by value.
            });
        }
    };

Finally, AMReX's expected OpenMP strategy for GPUs is to utilize OpenMP in
CPU regions to maintain multi-threaded parallelism on work that cannot be
offloaded efficiently, while using CUDA independently in GPU regions:
(MPI+OpenMP)+(MPI+CUDA).  This means OpenMP pragmas need to be maintained
when ``USE_CUDA=FALSE`` and turned off in locations CUDA is implemented
when ``USE_CUDA=TRUE``.

This can currently be implemented in preparation for an OpenMP strategy and
users are highly encouraged to do so now. This prevents having to track
down and label the appropriate OpenMP regions in the future and
clearly labels for readers that OpenMP and GPUs are not being used at the
same time.  OpenMP pragmas can be turned off using the conditional pragma
and :cpp:`Gpu::notInLaunchRegion()`, as shown below:

.. highlight:: c++

::

    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif

.. _sec:gpu:example:

An Example of Migrating to GPU
==============================

The nature of GPU programming poses difficulties for a number
of common AMReX patterns, such as the one below:

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

There are several issues in migrating this code to GPUs that need to
be addressed.  First, functions ``f1`` and ``f2`` have different
work regions (``vbx`` and ``gbx``, respectively) and there are data
dependencies between the two (``q``). This makes it difficult to put
them into a single GPU kernel, so two separate kernels will be
launched, one for each function.

As we have discussed in Section :ref:`sec:gpu:classes`, all
:cpp:`FArrayBox`\ es in the two :cpp:`MultiFab`\ s, ::cpp::`uin`
and :cpp:`uout` are in unified memory and avaiable on the GPUs.
But :cpp:`FArrayBox q` is in host memory.  Creating ``q`` as a 
managed object using the overloaded :cpp:`new` operator:

.. highlight:: c++

::

    FArrayBox* q = new FArrayBox;

does not solve the problem completely because GPU kernel calls are
asynchronous from CPU's point of view.  This creates a race
condition: GPU kernels in different iterations of :cpp:`MFIter`
will compete for access to ``q``.  One possible failure is a 
segfault when :cpp:`resize` changes the size of the ``q`` object
when the previous iteration of the loop is still using an old size.

Moving the line into the body of :cpp:`MFIter` loop will make ``q`` 
a variable local to each iteration, but it has a new issue.  When 
do we delete :cpp:`q`?  To the CPU, the resource of :cpp:`q` 
should be freed at the end of the scope, otherwise there will be 
a memory leak.  But at the end of the CPU scope, GPU kernels might 
still need it.

One way to fix this is put the temporary :cpp:`FArrayBox` objects in a
:cpp:`MultiFab` defined outside the loop.  This creates a separate
:cpp:`FArrayBox` for each loop iteration, eliminating the race
condition.  Another way is to use :cpp:`Gpu::AsyncFab` designed for 
this kind of situation.  The code below shows how :cpp:`Gpu::AsyncFab`
is used and how this MFIter loop can be rewritten for GPUs. 

.. highlight:: c++

::

   for (MFIter mfi(uin); mfi.isValid(); ++mfi)
   {
       const Box& vbx = mfi.validbox();              // f2 work domain
       const Box& gbx = amrex::grow(vbx,1);          // f1 work domain
       Gpu::AsyncFab q(gbx);                         // Local, GPU managed FArrayBox
       FArrayBox const* uinfab  = uin.fabPtr();      // Managed GPU capturable
       FArrayBox      * uoutfab = uout.fabPtr();     //   pointers to Multifab's FABs.

       AMREX_LAUNCH_DEVICE_LAMBDA ( gbx, tbx,        // f1 GPU launch 
       {
           f1(tbx, q.fab(), *uinfab);
       };

       AMrEX_LAMBDA_DEVICE_LAMBDA ( vbx, tbx,        // f2 GPU launch 
       {
           f2(tbx, *uoutfab, q.fab());
       });
   }                                                 // Implicit GPU barrier after
                                                     //   all iters are launched.

.. ===================================================================

.. _sec:gpu:assertion:

Assertions, Error Checking and Synchronization
================================================

To help debugging, we often use :cpp:`amrex::Assert` and
:cpp:`amrex::Abort`.  These functions are GPU safe and can be used in
GPU kernels.  However, implementing these functions requires additional
GPU registers, which will reduce overall performance.  Therefore, it
is preferred to implement such calls in debug mode only by wraping the
calls using ``#ifdef AMREX_DEBUG``. 

In CPU code, :cpp:`AMREX_GPU_ERROR_CHECK()` can be called
to check the health of previous GPU launches.  This call
looks up the return message from the most recently completed GPU
launch and aborts if it was not successful. Many kernel
launch macros as well as the :cpp:`MFIter` destructor include a call 
to :cpp:`AMREX_CPU_ERROR_CHECK()`. This prevents additional launches
from being called if a previous launch caused an error and ensures
all GPU launches within an :cpp:`MFIter` loop completed successfully
before continuing work.

However, due to asynchronicity, determining the source of the error 
can be difficult.  Even if GPU kernels launched earlier in the code 
result in a CUDA error, the error may not be output at a nearby call to
:cpp:`AMREX_GPU_ERROR_CHECK()` by the CPU.  When tracking down a CUDA
launch error, :cpp:`Gpu::Device::synchronize()` and 
:cpp:`Gpu::Device::streamSynchronize()` can be used to synchronize
the device or the CUDA stream, respectively, and track down the specific
launch that causes the error.

.. ===================================================================

Particle Support
================

.. _sec:gpu:particle:

AMReX's GPU particle support relies on Thrust, a parallel algorithms library maintained by
Nvidia. Thrust provides a GPU-capable vector container that is otherwise similar to the one
in the C++ Standard Template Library, along with associated sorting, searching, and prefix
summing operations. Combined with Cuda's unified memory, Thrust forms the basis of AMReX's
GPU support for particles. 

When compiled with ``USE_CUDA=TRUE``, AMReX places all its particle data in instances of
``thrust::device_vector`` that have been configured using a custom memory allocator using 
``cudaMallocManaged``. This means that the :cpp:`dataPtr` associated with particle data 
is managed and can be passed into GPU kernels, similar to the way it would be passed into
a Fortran subroutine in typical AMReX CPU code. As with the mesh data, these kernels can
be launched with a variety of approaches, including Cuda C / Fortran and OpenACC. An example
Fortran particle subroutine offloaded via OpenACC might look like the following:

.. highlight:: fortran

::

   subroutine push_position_boris(np, structs, uxp, uyp, uzp, gaminv, dt)

   use em_particle_module, only : particle_t
   use amrex_fort_module, only : amrex_real
   implicit none
   
   integer,          intent(in), value  :: np
   type(particle_t), intent(inout)      :: structs(np)
   real(amrex_real), intent(in)         :: uxp(np), uyp(np), uzp(np), gaminv(np)
   real(amrex_real), intent(in), value  :: dt
      
   integer                              :: ip

   !$acc parallel deviceptr(structs, uxp, uyp, uzp, gaminv)
   !$acc loop gang vector
   do ip = 1, np
       structs(ip)%pos(1) = structs(ip)%pos(1) + uxp(ip)*gaminv(ip)*dt
       structs(ip)%pos(2) = structs(ip)%pos(2) + uxp(ip)*gaminv(ip)*dt
       structs(ip)%pos(3) = structs(ip)%pos(3) + uxp(ip)*gaminv(ip)*dt
   end do
   !$acc end loop
   !$acc end parallel

   end subroutine push_position_boris
      
Note the use of the :fortran:`!$acc parallel deviceptr` clause to specify which data has been placed
in managed memory. This instructs OpenACC to treat those variables as if they already live on
the device, bypassing the usual copies. For a complete example of a particle code that has been ported
to GPUs using OpenACC, please see :cpp:`Tutorials/Particles/ElectromagneticPIC`. 
      
For portability, we have provided a set of Vector classes that wrap around the Thrust and
STL vectors. When ``USE_CUDA = FALSE``, these classes reduce to the normal :cpp:`amrex::Vector`.
When ``USE_CUDA = TRUE``, they have different meanings. :cpp:`Cuda::HostVector` is a wrapper
around :cpp:`thrust::host_vector`. :cpp:`Cuda::DeviceVector` is a wrapper around :cpp:`thrust::device_vector`,
while :cpp:`Cuda::ManagedDeviceVector` is a :cpp:`thrust::device_vector` that lives in managed memory.
These classes are useful when there are certain stages of an algorithm that will always
execute on either the host or the device. For example, the following code generates particles on
the CPU and copies them over to the GPU in one batch per tile:

.. highlight:: cpp

::

       for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
       {
           const Box& tile_box  = mfi.tilebox();      
           Cuda::HostVector<ParticleType> host_particles;
                           
           for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
           {
               < generate some particles... >
           }

           auto& particles = GetParticles(lev);
           auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
           auto old_size = particle_tile.GetArrayOfStructs().size();
           auto new_size = old_size + host_particles.size();
           particle_tile.resize(new_size);
           
           Cuda::thrust_copy(host_particles.begin(),
                             host_particles.end(),
                             particle_tile.GetArrayOfStructs().begin() + old_size);
        }

The following example shows how to use :cpp:`Cuda::DeviceVector`. Specifically, this code creates
temporary device vectors for the particle x, y, and z positions, and then copies from an Array-of-Structs
to a Struct-of-Arrays representation, all without copying any particle data off the GPU:

.. highlight:: cpp

::
   
   Cuda::DeviceVector<Real> xp, yp, zp;

   for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
   {
       pti.GetPosition(xp, yp, zp);

       < use xp, yp, zp... >
   }
           
Finally, AMReX's :cpp:`Redistribute()`, which moves particles back to the proper grids after their positions
have changed, has been ported to work on the GPU as well. It cannot be called from device code,
but it can be called on particles that reside on the device and it won't trigger any unified
memory traffic. As with :cpp:`MultiFab` data, the MPI portion of the particle redistribute is set
up to take advantange of the Cuda-aware MPI implementations available on platforms such as
ORNL's Summit and Summit-dev.

Performance Tips
================

.. _sec:gpu:performance:

Here are some helpful performance tips to keep in mind when working with AMReX for GPUs:

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
