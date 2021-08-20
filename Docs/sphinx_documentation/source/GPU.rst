.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:gpu:overview:

Overview of AMReX GPU Strategy
==============================

AMReX's GPU strategy focuses on providing performant GPU support
with minimal changes and maximum flexibility.  This allows
application teams to get running on GPUs quickly while allowing
long term performance tuning and programming model selection.  AMReX
uses the native programming language for GPUs: CUDA for NVIDIA, HIP
for AMD and DPC++ for Intel. This will be designated with ``CUDA/HIP/DPC++``
throughout the documentation.  However, application teams can also use
OpenACC or OpenMP in their individual codes.

At this time, AMReX does not support cross-native language compliation
(HIP for non-AMD systems and DPC++ for non Intel systems).  It may work with
a given version, but AMReX does not track or guarantee such functionality.

When running AMReX on a CPU system, the parallelization strategy is a
combination of MPI and OpenMP using tiling, as detailed in
:ref:`sec:basics:mfiter:tiling`. However, tiling is ineffective on GPUs
due to the overhead associated with kernel launching.  Instead,
efficient use of the GPU's resources is the primary concern.  Improving
resource efficiency allows a larger percentage of GPU threads to work
simultaneously, increasing effective parallelism and decreasing the time
to solution.

When running on CPUs, AMReX uses an ``MPI+X`` strategy where the ``X``
threads are used to perform parallelization techniques, like tiling.
The most common ``X`` is ``OpenMP``.  On GPUs, AMReX requires ``CUDA/HIP/DPC++``
and can be further combined with other parallel GPU languages, including
``OpenACC`` and ``OpenMP``, to control the offloading of subroutines
to the GPU.  This ``MPI+CUDA+X`` GPU strategy has been developed
to give users the maximum flexibility to find the best combination of
portability, readability and performance for their applications.

Presented here is an overview of important features of AMReX's GPU strategy.
Additional information that is required for creating GPU applications is
detailed throughout the rest of this chapter:

- Each MPI rank offloads its work to a single GPU. ``(MPI ranks == Number of GPUs)``

- Calculations that can be offloaded efficiently to GPUs use GPU threads
  to parallelize over a valid box at a time.  This is done by launching over
  a large number GPU threads that only work on a few cells each. This work
  distribution is illustrated in :numref:`fig:gpu:threads`.

.. |a| image:: ./GPU/gpu_2.png
       :width: 100%

.. |b| image:: ./GPU/gpu_3.png
       :width: 100%

.. _fig:gpu:threads:

.. table:: Comparison of OpenMP and GPU work distribution. Pictures provided by Mike Zingale and the CASTRO team.

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a|                          |                        |b|                           |
   +-----------------------------------------------------+------------------------------------------------------+
   |   OpenMP tiled box.                                 |   GPU threaded box.                                  |
   |   OpenMP threads break down the valid box           |   Each GPU thread works on a few cells of the        |
   |   into two large boxes (blue and orange).           |   valid box. This example uses one cell per          |
   |   The lo and hi of one tiled box are marked.        |   thread, each thread using a box with lo = hi.      |
   +-----------------------------------------------------+------------------------------------------------------+

- C++ macros and GPU extended lambdas are used to provide performance
  portability while making the code as understandable as possible to
  science-focused code teams.

- AMReX utilizes GPU managed memory to automatically handle memory
  movement for mesh and particle data.  Simple data structures, such
  as :cpp:`IntVect`\s can be passed by value and complex data structures, such as
  :cpp:`FArrayBox`\es, have specialized AMReX classes to handle the
  data movement for the user.  Tests have shown CUDA managed memory
  to be efficient and reliable, especially when applications remove
  any unnecessary data accesses.

- Application teams should strive to keep mesh and particle data structures
  on the GPU for as long as possible, minimizing movement back to the CPU.
  This strategy lends itself to AMReX applications readily; the mesh and
  particle data can stay on the GPU for most subroutines except for
  of redistribution, communication and I/O operations.

- AMReX's GPU strategy is focused on launching GPU kernels inside AMReX's
  :cpp:`MFIter` and :cpp:`ParIter` loops.  By performing GPU work within
  :cpp:`MFIter` and :cpp:`ParIter` loops, GPU work is isolated to independent
  data sets on well-established AMReX data objects, providing consistency and safety
  that also matches AMReX's coding methodology.  Similar tools are also available for
  launching work outside of AMReX loops.

- AMReX further parallelizes GPU applications by utilizing streams.
  Streams guarantee execution order of kernels within the same stream, while
  allowing different streams to run simultaneously. AMReX places each iteration
  of :cpp:`MFIter` loops on separate streams, allowing each independent
  iteration to be run simultaneously and sequentially, while maximizing GPU usage.

  The AMReX implementation of streams is illustrated in :numref:`fig:gpu:streams`.
  The CPU runs the first iteration of the MFIter loop (blue), which contains three
  GPU kernels.  The kernels begin immediately in GPU Stream 1 and run in the same
  order they were added. The second (red) and third (green) iterations are similarly
  launched in Streams 2 and 3. The fourth (orange) and fifth (purple) iterations
  require more GPU resources than remain, so they have to wait until resources are
  freed before beginning. Meanwhile, after all the loop iterations are launched, the
  CPU reaches a synchronize in the MFIter's destructor and waits for all GPU launches
  to complete before continuing.

- The Fortran interface of AMReX does not currently have GPU support.  AMReX recommends
  porting Fortran code to C++ when coding for GPUs.

.. raw:: latex

   \begin{center}

.. _fig:gpu:streams:

.. figure:: ./GPU/Streams.png

   Timeline illustration of GPU streams. Illustrates the case of an
   MFIter loop of five iterations with three GPU kernels each being
   ran with three GPU streams.

.. raw:: latex

   \end{center}

.. _sec:gpu:build:

Building GPU Support
====================

Building with GNU Make
----------------------

To build AMReX with GPU support, add ``USE_CUDA=TRUE``, ``USE_HIP=TRUE`` or
``USE_DPCPP=TRUE`` to the ``GNUmakefile`` or as a command line argument.

AMReX does not require OpenACC, but application codes
can use them if they are supported by the compiler.  For OpenACC support, add
``USE_ACC=TRUE``.  PGI, Cray and GNU compilers support OpenACC.  Thus,
for OpenACC, you must use ``COMP=pgi``, ``COMP=cray`` or ``COMP=gnu``.

Currently, only IBM is supported with OpenMP offloading. To use OpenMP
offloading, make with ``USE_OMP_OFFLOAD=TRUE``.

Compiling AMReX with CUDA requires compiling the code through NVIDIA's
CUDA compiler driver in addition to the standard compiler.  This driver
is called ``nvcc`` and it requires a host compiler to work through.
The default host compiler for NVCC is GCC even if ``COMP`` is set to
a different compiler.  One can change this by setting ``NVCC_HOST_COMP``.
For example, ``COMP=pgi`` alone will compile C/C++ codes with NVCC/GCC
and Fortran codes with PGI, and link with PGI.  Using ``COMP=pgi`` and
``NVCC_HOST_COMP=pgi`` will compile C/C++ codes with PGI and NVCC/PGI.

You can use ``Tutorials/Basic/HelloWorld_C`` to test your programming
environment.  For example, building with:

.. highlight:: console

::

   make COMP=gnu USE_CUDA=TRUE

should produce an executable named ``main3d.gnu.DEBUG.CUDA.ex``.  You
can run it and that will generate results like:

.. highlight:: console

::

   $ ./main3d.gnu.DEBUG.CUDA.ex
   Initializing CUDA...
   CUDA initialized with 1 GPU
   AMReX (19.06-404-g0455b168b69c-dirty) initialized
   Hello world from AMReX version 19.06-404-g0455b168b69c-dirty
   Total GPU global memory (MB): 6069
   Free  GPU global memory (MB): 5896
   [The         Arena] space (MB): 4552
   [The Managed Arena] space (MB): 8
   [The  Pinned Arena] space (MB): 8
   AMReX (19.06-404-g0455b168b69c-dirty) finalized

Building with CMake
-------------------

Enabling CUDA support
^^^^^^^^^^^^^^^^^^^^^

To build AMReX with CUDA support in CMake, add ``-DAMReX_GPU_BACKEND=CUDA`` to the
``cmake`` invocation. For a full list of CUDA-specific configuration options,
check the :ref:`table <tab:cmakecudavar>` below.

.. raw:: latex

   \begin{center}

.. _tab:cmakecudavar:

.. table:: AMReX CUDA-specific build options

   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | Variable Name                | Description                                     | Default     | Possible values |
   +==============================+=================================================+=============+=================+
   | AMReX_CUDA_ARCH              |  CUDA target architecture                       | Auto        | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_FASTMATH          |  Enable CUDA fastmath library                   | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_BACKTRACE         |  Host function symbol names (e.g. cuda-memcheck)| Auto        | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_COMPILATION_TIMER |  CSV table with time for each compilation phase | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_DEBUG             |  Device debug information (optimizations: off)  | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_ERROR_CAPTURE_THIS|  Error if a CUDA lambda captures a class' this  | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_ERROR_CROSS       |  Error if a host function is called from a host | NO          | YES, NO         |
   |  _EXECUTION_SPACE_CALL       |   device function                               |             |                 |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_KEEP_FILES        |  Keep intermediately files (folder: nvcc_tmp)   | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_LTO               |  Enable CUDA link-time-optimization             | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_MAX_THREADS       |  Max number of CUDA threads per block           | 256         | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_MAXREGCOUNT       |  Limits the number of CUDA registers available  | 255         | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_PTX_VERBOSE       |  Verbose code generation statistics in ptxas    | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_SHOW_CODELINES    |  Source information in PTX (optimizations: on)  | Auto        | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_SHOW_LINENUMBERS  |  Line-number information (optimizations: on)    | Auto        | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_CUDA_WARN_CAPTURE_THIS |  Warn if a CUDA lambda captures a class' this   | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
.. raw:: latex

   \end{center}


The target architecture to build for can be specified via the configuration option
``-DAMReX_CUDA_ARCH=<target-architecture>``, where ``<target-architecture>`` can be either
the name of the NVIDIA GPU generation, i.e. ``Turing``, ``Volta``, ``Ampere``, ``...`` , or its
`compute capability <https://developer.nvidia.com/cuda-gpus>`_, i.e. ``10.0``, ``9.0``,  ``...`` .
For example, on Cori GPUs you can specify the architecture as follows:

.. highlight:: console

::

   cmake [options] -DAMReX_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=Volta /path/to/amrex/source


If no architecture is specified, CMake will default to the architecture defined in the
*environment variable* ``AMREX_CUDA_ARCH`` (note: all caps).
If the latter is not defined, CMake will try to determine which GPU architecture is supported by the system.
If more than one is found, CMake will build for all of them.
If autodetection fails, a list of "common" architectures is assumed.
`Multiple CUDA architectures <https://cmake.org/cmake/help/latest/module/FindCUDA.html#commands>`__ can also be set manually as semicolon-separated list, e.g. ``-DAMReX_CUDA_ARCH=7.0;8.0``.
Building for multiple CUDA architectures will generally result in a larger library and longer build times.

**Note that AMReX supports NVIDIA GPU architectures with compute capability 6.0 or higher and
CUDA Toolkit version 9.0 or higher.**

In order to import the CUDA-enabled AMReX library into your CMake project, you need to include
the following code into the appropriate CMakeLists.txt file:

.. highlight:: console

::

   # Find CUDA-enabled AMReX installation
   find_package(AMReX REQUIRED CUDA)


If instead of using an external installation of AMReX you prefer to include AMReX as a subproject
in your CMake setup, we strongly encourage you to use the ``AMReX_SetupCUDA`` module as shown below
if the CMake version is less than 3.20:

.. highlight:: console

::

   # Enable CUDA in your CMake project
   enable_language(CUDA)

   # Include the AMReX-provided CUDA setup module -- OBSOLETE with CMake >= 3.20
   if(CMAKE_VERSION VERSION_LESS 3.20)
       include(AMReX_SetupCUDA)
   endif()

   # Include AMReX source directory ONLY AFTER the two steps above
   add_subdirectory(/path/to/amrex/source/dir)



To ensure consistency between CUDA-enabled AMReX and any CMake target that links against it,
we provide the helper function ``setup_target_for_cuda_compilation()``:


.. highlight:: console

::

   # Set all sources for my_target
   target_sources(my_target source1 source2 source3 ...)

   # Setup my_target to be compiled with CUDA and be linked against CUDA-enabled AMReX
   # MUST be done AFTER all sources have been assigned to my_target
   setup_target_for_cuda_compilation(my_target)

   # Link against amrex
   target_link_libraries(my_target AMReX::amrex)



Enabling HIP support (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build AMReX with HIP support in CMake, add
``-DAMReX_GPU_BACKEND=HIP -DAMReX_AMD_ARCH=<target-arch> -DCMAKE_CXX_COMPILER=<your-hip-compiler>``
to the ``cmake`` invocation.
If you don't need Fortran features (``AMReX_FORTRAN=OFF``), it is recommended to use AMD's ``clang++`` as the HIP compiler.
(Please see these issues for reference in rocm/HIP <= 4.2.0
`[1] <https://github.com/ROCm-Developer-Tools/HIP/issues/2275>`__
`[2] <https://github.com/AMReX-Codes/amrex/pull/2031>`__.)

In AMReX CMake, the HIP compiler is treated as a special C++ compiler and therefore
the standard CMake variables used to customize the compilation process for C++,
for example ``CMAKE_CXX_FLAGS``, can be used for HIP as well.


Since CMake does not support autodetection of HIP compilers/target architectures
yet, ``CMAKE_CXX_COMPILER`` must be set to a valid HIP compiler, i.e. ``clang++`` or ``hipcc`` or ``nvcc``,
and ``AMReX_AMD_ARCH`` to the target architecture you are building for.
Thus **AMReX_AMD_ARCH and CMAKE_CXX_COMPILER are required user-inputs when AMReX_GPU_BACKEND=HIP**.
We again read also an *environment variable*: ``AMREX_AMD_ARCH`` (note: all caps) and the C++ compiler can be hinted as always, e.g. with ``export CXX=$(which clang++)``.
Below is an example configuration for HIP on Tulip:

.. highlight:: console

::

   cmake -S . -B build -DAMReX_GPU_BACKEND=HIP -DCMAKE_CXX_COMPILER=$(which clang++) -DAMReX_AMD_ARCH="gfx906;gfx908"  # [other options]
   cmake --build build -j 6


Enabling SYCL support (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build AMReX with SYCL support in CMake, add
``-DAMReX_GPU_BACKEND=SYCL -DCMAKE_CXX_COMPILER=<your-sycl-compiler>``
to the ``cmake`` invocation.
For a full list of SYCL-specific configuration options,
check the :ref:`table <tab:cmakesyclvar>` below.


In AMReX CMake, the SYCL compiler is treated as a special C++ compiler and therefore
the standard CMake variables used to customize the compilation process for C++,
for example ``CMAKE_CXX_FLAGS``, can be used for DPCPP as well.


Since CMake does not support autodetection of SYCL compilers yet,
``CMAKE_CXX_COMPILER`` must be set to a valid SYCL compiler. i.e. ``dpcpp``.
Thus **CMAKE_CXX_COMPILER is a required user-input when AMReX_GPU_BACKEND=SYCL**.
At this time, **the only supported SYCL compiler is dpcpp**.
Below is an example configuration for SYCL:

.. highlight:: console

::

   cmake -DAMReX_GPU_BACKEND=SYCL -DCMAKE_CXX_COMPILER=$(which dpcpp)  [other options] /path/to/amrex/source


.. raw:: latex

   \begin{center}

.. _tab:cmakesyclvar:

.. table:: AMReX SYCL-specific build options

   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | Variable Name                | Description                                     | Default     | Possible values |
   +==============================+=================================================+=============+=================+
   | AMReX_DPCPP_AOT              | Enable DPCPP ahead-of-time compilation          | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_DPCPP_SPLIT_KERNEL     | Enable DPCPP kernel splitting                   | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | AMReX_DPCPP_ONEDPL           | Enable DPCPP's oneDPL algorithms                | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
.. raw:: latex

   \end{center}



.. ===================================================================

.. _sec:gpu:namespace:

Gpu Namespace and Macros
========================

Most GPU related classes and functions are in ``namespace Gpu``,
which is inside ``namespace amrex``. For example, the GPU configuration
class ``Device`` can be referenced to at ``amrex::Gpu::Device``.

For portability, AMReX defines some macros for CUDA function qualifiers
and they should be preferred to allow execution with ``USE_CUDA=FALSE``.
These include:

.. highlight:: c++

::

   #define AMREX_GPU_HOST        __host__
   #define AMREX_GPU_DEVICE      __device__
   #define AMREX_GPU_GLOBAL      __global__
   #define AMREX_GPU_HOST_DEVICE __host__ __device__

Note that when AMReX is not built with ``CUDA/HIP/DPC++``,
these macros expand to empty space.

When AMReX is compiled with ``USE_CUDA=TRUE``, the preprocessor
macros ``AMREX_USE_CUDA`` and ``AMREX_USE_GPU`` are defined for
conditional programming.  When AMReX is compiled with
``USE_ACC=TRUE``, ``AMREX_USE_ACC`` is defined.  When AMReX is
compiled with ``USE_OMP_OFFLOAD=TRUE``, ``AMREX_USE_OMP_OFFLOAD`` is
defined.

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
CUDA, all :cpp:`Arena`\ s use standard :cpp:`new` and :cpp:`delete`
operators. With CUDA, the :cpp:`Arena`\ s each allocate with a
specific type of GPU memory:

.. raw:: latex

    \begin{center}

.. _tab:gpu:arena:

.. table:: Memory Arenas

    +---------------------+----------------------------+
    | Arena               |        Memory Type         |
    +=====================+============================+
    | The_Arena()         |  managed or device memory  |
    +---------------------+----------------------------+
    | The_Managed_Arena() |  managed memory            |
    +---------------------+----------------------------+
    | The_Pinned_Arena()  |  pinned memory             |
    +---------------------+----------------------------+

.. raw:: latex

    \end{center}

The Arena object returned by these calls provides access
to two functions:

.. highlight:: c++

::

   void* alloc (std::size_t sz);
   void free (void* p);

:cpp:`The_Arena()` is used for memory allocation of data in
:cpp:`BaseFab`.  By default, it allocates managed memory.  This can be changed with
a boolean runtime parameter ``amrex.the_arena_is_managed``.
Therefore the data in a :cpp:`MultiFab` is placed in
managed memory by default and is accessible from both CPU host and GPU device.
This allows application codes to develop their GPU capability
gradually. The behavior of :cpp:`The_Managed_Arena()` likewise depends on the
``amrex.the_arena_is_managed`` parameter. If ``amrex.the_arena_is_managed=0``,
:cpp:`The_Managed_Arena()` is a separate pool of managed memory. If
``amrex.the_arena_is_managed=1``, :cpp:`The_Managed_Arena()` is simply aliased
to :cpp:`The_Arena()` to reduce memory fragmentation.

If you want to print out the current memory usage
of the Arenas, you can call :cpp:`amrex::Arena::PrintUsage()`.
When AMReX is built with SUNDIALS turned on, :cpp:`amrex::sundials::The_SUNMemory_Helper()`
can be provided to SUNDIALS data structures so that they use the appropriate
Arena object when allocating memory. For example, it can be provided to the
SUNDIALS CUDA vector:

.. highlight:: c++

::

  N_Vector x = N_VNewWithMemHelp_Cuda(size, use_managed_memory, *The_SUNMemory_Helper());


.. ===================================================================

.. _sec:gpu:classes:

GPU Safe Classes and Functions
==============================

AMReX GPU work takes place inside of MFIter and particle loops.
Therefore, there are two ways classes and functions have been modified
to interact with the GPU:

1. A number of functions used within these loops are labelled using
``AMREX_GPU_HOST_DEVICE`` and can be called on the device. This includes member
functions, such as :cpp:`IntVect::type()`, as well as non-member functions,
such as :cpp:`amrex::min` and :cpp:`amrex::max`. In specialized cases,
classes are labeled such that the object can be constructed, destructed
and its functions can be implemented on the device, including ``IntVect``.

2. Functions that contain MFIter or particle loops have been rewritten
to contain device launches. For example, the :cpp:`FillBoundary`
function cannot be called from device code, but calling it from
CPU will launch GPU kernels if AMReX is compiled with GPU support.

Necessary and convenient AMReX functions and objects have been given a device
version and/or device access.

In this section, we discuss some examples of AMReX device classes and functions
that are important for programming GPUs.


GpuArray, Array1D, Array2D, and Array3D
---------------------------------------

As we have mentioned in :ref:`sec:basics:vecandarr`, :cpp:`std::array`
cannot be used in device code, whereas :cpp:`GpuArray`,
:cpp:`Array1D`, :cpp:`Array2D`, and :cpp:`Array3D` are trivial types
that work on both host and device. They can be used whenever a fixed size array
needs to be passed to the GPU or created on GPU.  A variety of
functions in AMReX return :cpp:`GpuArray` and they can be
lambda-captured to GPU code. For example,
:cpp:`GeometryData::CellSizeArray()`, :cpp:`GeometryData::InvCellSizeArray()`
and :cpp:`Box::length3d()` all return :cpp:`GpuArray`\s.

.. _sec:gpu:classes:asyncarray:


AsyncArray
----------

Where the :cpp:`GpuArray` is a statically-sized array designed to be
passed by value onto the device, :cpp:`AsyncArray` is a
dynamically-sized array container designed to work between the CPU and
GPU. :cpp:`AsyncArray` stores a CPU pointer and a GPU pointer and
coordinates the movement of an array of objects between the two.  It
can take initial values from the host and move them to the device.  It
can copy the data from device back to host.  It can also be used as
scratch space on device.

The call to delete the memory is added to the GPU stream as a callback
function in the destructor of :cpp:`AsyncArray`. This guarantees the
memory allocated in :cpp:`AsyncArray` continues to exist after the
:cpp:`AsyncArray` object is deleted when going out of scope until
after all GPU kernels in the stream are completed without forcing the
code to synchronize. The resulting :cpp:`AsyncArray` class is
"async-safe", meaning it can be safely used in asynchronous code
regions that contain both CPU work and GPU launches, including
:cpp:`MFIter` loops.

:cpp:`AsyncArray` is also portable. When built without ``USE_CUDA``, the
object only stores and handles the CPU version of the data.

An example using :cpp:`AsyncArray` is given below,

.. highlight:: c++

::

    Real h_s = 0.0;
    AsyncArray<Real> aa_s(&h_s, 1);  // Build AsyncArray of size 1
    Real* d_s = aa_s.data();         // Get associated device pointer

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Vector<Real> h_v = a_cpu_function();
        AsyncArray<Real> aa_v1(h_v.data(), h_v.size());
        Real* d_v1 = aa_v1.data();  // A device copy of the data

        std::size_t n = ...;
        AsyncArray<Real> aa_v2(n);  // Allocate temporary space on device
        Real* d_v2 = aa_v2.data();  // A device pointer to uninitialized data

        ... // gpu kernels using the data pointed by d_v1 and atomically
            // updating the data pointed by d_s.
            // d_v2 can be used as scratch space and for pass data
            // between kernels.

        // If needed, we can copy the data back to host using
        // AsyncArray::copyToHost(host_pointer, number_of_elements);

        // At the end of each loop the compiler inserts a call to the
        // destructor of aa_v* on cpu.  Objects aa_v* are deleted, but
        // their associated memory pointed by d_v* is not deleted
        // immediately until the gpu kernels in this loop finish.
    }

    aa_s.copyToHost(&h_s, 1); // Copy the value back to host

Gpu Vectors
-----------

AMReX also provides a number of dynamic vectors for use with GPU kernels.
These are configured to use the different AMReX memory Arenas, as
summarized below. By using the memory Arenas, we can avoid expensive
allocations and deallocations when (for example) resizing vectors.

.. raw:: latex

    \begin{center}

.. _tab:gpu:gpuvectors:

.. table:: Memory Arenas Associated with each Gpu Vector

    +----------------+----------------------+
    | Vector         | Arena                |
    +================+======================+
    | DeviceVector   | The_Arena()          |
    +----------------+----------------------+
    | HostVector     | The_Pinned_Arena()   |
    +----------------+----------------------+
    | ManagedVector  | The_Managed_Arena()  |
    +----------------+----------------------+

.. raw:: latex

    \end{center}

These classes behave identically to an
:cpp:`amrex::Vector`, (see :ref:`sec:basics:vecandarr`), except that they
can only hold "plain-old-data" objects (e.g. Reals, integers, amrex Particles,
etc... ). If you want a resizable vector that doesn't use a memory Arena,
simply use :cpp:`amrex::Vector`.

Note that, even if the data in the vector is  managed and available on GPUs,
the member functions of e.g. :cpp:`Gpu::ManagedVector` are not.
To use the data on the GPU, it is necessary to pass the underlying data pointer
in to the GPU kernels. The managed data pointer can be accessed using the :cpp:`data()`
member function.

Be aware: resizing of dynamically allocated memory on the GPU is unsupported.
All resizing of the vector should be done on the CPU, in a manner that avoids
race conditions with concurrent GPU kernels.

Also note: :cpp:`Gpu::ManagedVector` is not async-safe.  It cannot be safely
constructed inside of an MFIter loop with GPU kernels and great care should
be used when accessing :cpp:`Gpu::ManagedVector` data on GPUs to avoid race
conditions.

amrex::min and amrex::max
-------------------------

GPU versions of ``std::min`` and ``std::max`` are not provided in CUDA.
So, AMReX provides a templated :cpp:`min` and :cpp:`max` with host and
device versions to allow functionality on GPUs. Invoke the explicitly
namespaced :cpp:`amrex::min(A, B)` or :cpp:`amrex::max(x, y)` to use the
GPU safe implementations. These functions are variadic, so they can take
any number of arguments and can be invoked with any standard data type.


MultiFab Reductions
-------------------

AMReX provides functions for performing standard reduction operations on
:cpp:`MultiFabs`, including :cpp:`MultiFab::sum` and :cpp:`MultiFab::max`.
When ``USE_CUDA=TRUE``, these functions automatically implement the
corresponding reductions on GPUs in an efficient manner.

Function templates :cpp:`amrex::ReduceSum`, :cpp:`amrex::ReduceMin` and
:cpp:`amrex::ReduceMax` can be used to implement user-defined reduction
functions over :cpp:`MultiFab`\ s. These same templates are implemented
in the :cpp:`MultiFab` functions, so they can be used as a reference to
build a custom reduction. For example, the :cpp:`MultiFab:Dot`
implementation is reproduced here:

.. highlight:: c++

::

    Real MultiFab::Dot (const MultiFab& x, int xcomp,
                        const MultiFab& y, int ycomp,
                        int numcomp, int nghost, bool local) {
        Real sm = amrex::ReduceSum(x, y, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        {
            return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp);
        });

        if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

        return sm;
    }

:cpp:`amrex::ReduceSum` takes two :cpp:`MultiFab`\ s, ``x`` and ``y`` and
returns the sum of the value returned from the given lambda function.
In this case, :cpp:`BaseFab::dot` is returned, yielding a sum of the
dot product of each local pair of :cpp:`BaseFab`\ s. Finally,
:cpp:`ParallelAllReduce` is used to sum the dot products across all
MPI ranks and return the total dot product of the two
:cpp:`MultiFab`\ s.

To implement a different reduction, replace the code block inside the
lambda function with the operation that should be applied, being sure
to return the value to be summed, minimized, or maximized.  The reduction
templates have a few different interfaces to accommodate a variety of
reductions.  The :cpp:`amrex::ReduceSum` reduction template has varieties
that take either one, two or three ::cpp:`MultiFab`\ s.
:cpp:`amrex::ReduceMin` and :cpp:`amrex::ReduceMax` can take either one
or two.


Box, IntVect and IndexType
--------------------------

In AMReX, :cpp:`Box`, :cpp:`IntVect` and :cpp:`IndexType`
are classes for representing indices.  These classes and most of
their member functions, including constructors and destructors,
have both host and device versions.  They can be used freely
in device code.


Geometry
--------

AMReX's :cpp:`Geometry` class is not a GPU safe class.  However, we often need
to use geometric information such as cell size and physical coordinates
in GPU kernels.  We can use the following member functions and pass
the returned values to GPU kernels:

.. highlight:: c++

::

    GpuArray<Real,AMREX_SPACEDIM> ProbLoArray () const noexcept;
    GpuArray<Real,AMREX_SPACEDIM> ProbHiArray () const noexcept;
    GpuArray<int,AMREX_SPACEDIM> isPeriodicArray () const noexcept;
    GpuArray<Real,AMREX_SPACEDIM> CellSizeArray () const noexcept;
    GpuArray<Real,AMREX_SPACEDIM> InvCellSizeArray () const noexcept;

Alternatively, we can copy the data into a GPU safe class that can be
passed by value to GPU kernels. This class is called
:cpp:`GeometryData`, which is created by calling
:cpp:`Geometry::data()`.  The accessor functions of
:cpp:`GeometryData` are identical to :cpp:`Geometry`.

.. _sec:gpu:classes:basefab:

BaseFab, FArrayBox, IArrayBox
-----------------------------

:cpp:`BaseFab<T>`, :cpp:`IArrayBox` and :cpp:`FArrayBox` have some GPU
support.  They cannot be constructed in device code unless they are
constructed as an alias to :cpp:`Array4`.  Many of their member
functions can be used in device code as long as they have been
constructed in device memory. Some of the device member functions
include :cpp:`array`, :cpp:`dataPtr`, :cpp:`box`, :cpp:`nComp`, and
:cpp:`setVal`.

All :cpp:`BaseFab<T>` objects in :cpp:`FabArray<FAB>` are allocated in
CPU memory, including :cpp:`IArrayBox` and :cpp:`FArrayBox`, which are
derived from :cpp:`BaseFab`, although the array data contained are
allocated in managed memory.  We cannot pass a :cpp:`BaseFab` object by
value because they do not have copy constructor.  However, we can make
an :cpp:`Array4` using member function :cpp:`BaseFab::array()`, and pass it
by value to GPU kernels. In GPU device code, we can use :cpp:`Array4`
or, if necessary, we can make an alias :cpp:`BaseFab` from an
:cpp:`Array4`.  For example,

.. highlight:: c++

::

    AMREX_GPU_HOST_DEVICE void g (FArrayBox& fab) { ... }

    AMREX_GPU_HOST_DEVICE void f (Box const& bx, Array4<Real> const& a)
    {
      FArrayBox fab(a,bx.ixType());
      g(fab);
    }

.. _sec:gpu:classes:elixir:

Elixir
------

We often have temporary :cpp:`FArrayBox`\ es in :cpp:`MFIter` loops.
These objects go out of scope at the end of each iteration.  Because
of the asynchronous nature of GPU kernel execution, their destructors
might get called before their data are used on GPU.  :cpp:`Elixir` can
be used to extend the life of the data.  For example,

.. highlight:: c++

::

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      FArrayBox tmp_fab(bx, numcomps);
      Elixir tmp_eli = tmp_fab.elixir();
      Array4<Real> const& tmp_arr = tmp_fab.array();

      // GPU kernels using the temporary
    }

Without :cpp:`Elixir`, the code above will likely cause memory errors
because the temporary :cpp:`FArrayBox` is deleted on cpu before the
gpu kernels use its memory.  With :cpp:`Elixir`, the ownership of the
memory is transferred to :cpp:`Elixir` that is guaranteed to be
async-safe.

Async Arena
-----------

CUDA 11.2 has introduced a new feature, stream-ordered CUDA memory
allocator.  This feature enables AMReX to solve the temporary memory
allocation and deallocation issue discussed above using a memory pool.
Instead of using :cpp:`Elixir`, we can write code like below,

.. highlight:: c++

::

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      FArrayBox tmp_fab(bx, numcomps, The_Async_Arena());
      Array4<Real> const& tmp_arr = tmp_fab.array();
      FArrayBox tmp_fab_2;
      tmp_fab_2.resize(bx, numcomps, The_Async_Arena());

      // GPU kernels using the temporary
    }

This is now the recommended way because it's usually more efficient than
:cpp:`Elixir`.  Note that the code above works for CUDA older than 11.2, HIP
and DPC++ as well, and it's equivalent to using :cpp:`Elixir` in these
cases.  By default, the release threshold for the memory pool is unlimited.
One can adjust it with :cpp:`ParmParse` parameter,
``amrex.the_async_arena_release_threshold``.

.. _sec:gpu:launch:

Kernel Launch
=============

In this section, how to offload work to the GPU will be demonstrated.
AMReX supports offloading work with CUDA, OpenACC, or OpenMP.

When using CUDA, AMReX provides users with portable C++ function calls or
C++ macros that launch a user-defined lambda function.  When compiled without CUDA,
the lambda function is ran on the CPU. When compiled with CUDA, the launch function
prepares and launches the lambda function on the GPU. The preparation includes
calculating the appropriate number of blocks and threads, selecting the CUDA stream
and defining the appropriate work chunk for each CUDA thread.

When using OpenACC or OpenMP offloading pragmas, the users add the appropriate
pragmas to their work loops and functions to offload to the GPU.  These work
in conjunction with AMReX's internal CUDA-based memory management, described
earlier, to ensure the required data is available on the GPU when the offloaded
function is executed.

The available launch schema are presented here in three categories: launching
nested loops over Boxes or 1D arrays, launching generic work and launching using
OpenACC or  OpenMP pragmas. The latest versions of the examples used in this section
of the documentation can be found in the AMReX source code in the `Launch`_ tutorials.
Users should also refer to Chapter :ref:`Chap:Basics` as needed for information about basic
AMReX classes.

.. _`Launch`: https://amrex-codes.github.io/amrex/tutorials_html/GPU_Tutorial.html#launch


AMReX also recommends writing primary floating point operation kernels
in C++ using AMReX's :cpp:`Array4` object syntax.  It provides a
multi-dimensional array syntax, similar in appearance to Fortran,
while maintaining performance.  The details can be found in
:ref:`Array4 <sec:basics:array4>` and :ref:`C++ Kernel
<sec:basics:cppkernel>`.

.. Overview table???

.. _sec:gpu:for:

Launching C++ nested loops
--------------------------

The most common AMReX work construct is a set of nested loops over
the cells in a box. AMReX provides C++ functions and macro equivalents to port nested
loops efficiently onto the GPU.  There are 3 different nested loop GPU
launches: a 4D launch for work over a box and a number of components, a 3D
launch for work over a box and a 1D launch for work over a number of arbitrary elements.
Each of these launches provides a performance portable set of nested loops for
both CPU and GPU applications.

These loop launches should only be used when each iteration of the
nested loop is independent of other iterations.  Therefore, these
launches have been marked with ``AMREX_PRAGMA_SIMD`` when using the
CPU and they should only be used for ``simd``-capable nested loops.
Calculations that cannot vectorize should be rewritten wherever
possible to allow efficient utilization of GPU hardware.

However, it is important for applications to use these launches whenever appropriate
because they contain optimizations for both CPU and GPU variations of nested
loops.  For example, on the GPU the spatial coordinate loops are reduced to a single
loop and the component loop is moved to these inner most loop.  AMReX's launch functions
apply the appropriate optimizations for ``USE_CUDA=TRUE`` and ``USE_CUDA=FALSE`` in a
compact and readable format.

AMReX also provides a variation of the launch function that is implemented as a
C++ macro.  It behaves identically to the function, but hides the lambda function
from to the user.  There are some subtle differences between the two implementations,
that will be discussed.  It is up to the user to select which version they would like
to use.  For simplicity, the function variation will be discussed throughout the rest of
this documentation, however all code snippets will also include the macro variation
for reference.

A 4D example of the launch function, :cpp:`amrex::ParallelFor`, is given here:

.. highlight:: c++

::

    int ncomp = mf.nComp();
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& fab = mf.array(mfi);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            fab(i,j,k,n) += 1.;
        });

        /* MACRO VARIATION:
        /
        /   AMREX_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        /   {
        /       fab(i,j,k,n) += 1.;
        /   });
        */
    }

This code works whether it is compiled for GPUs or CPUs. :cpp:`TilingIfNotGPU()`
returns ``false`` in the GPU case to turn off tiling and maximize the amount of
work given to the GPU in each launch. When tiling is off, :cpp:`tilebox()`
returns the :cpp:`validbox()`.  The :cpp:`BaseFab::array()` function returns a
lightweight :cpp:`Array4` object that defines access to the underlying :cpp:`FArrayBox`
data.  The :cpp:`Array4`\s is then captured by the C++ lambda functions defined in the
launch function.

``amrex::ParallelFor()`` expands into different variations of a quadruply-nested
:cpp:`for` loop depending dimensionality and whether it is being implemented on CPU or GPU.
The best way to understand this macro is to take a look at the 4D :cpp:`amrex::ParallelFor`
that is implemented when ``USE_CUDA=FALSE``. A simplified version is reproduced here:

.. highlight:: c++

::

    void ParallelFor (Box const& box, int ncomp, /* LAMBDA FUNCTION */)
    {
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        for (int n = 0; n < ncomp; ++n) {
            for (int z = lo.z; z <= hi.z; ++z) {
            for (int y = lo.y; y <= hi.y; ++y) {
            AMREX_PRAGMA_SIMD
            for (int x = lo.x; x <= hi.x; ++x) {
                /* LAUNCH LAMBDA FUNCTION (x,y,z,n) */
            }}}
        }
    }

:cpp:`amrex::ParallelFor` takes a :cpp:`Box` and a number of components, which define the bounds
of the quadruply-nested :cpp:`for` loop, and a lambda function to run on each iteration of the
nested loop.  The lambda function takes the loop iterators as parameters, allowing the current
cell to be indexed in the lambda.  In addition to the loop indices, the lambda function captures
any necessary objects defined in the local scope.

CUDA lambda functions can only capture by value, as the information
must be able to be copied onto the device.  In this example, the
lambda function captures a :cpp:`Array4` object, ``fab``, that defines
how to access the :cpp:`FArrayBox`.  The macro uses ``fab`` to
increment the value of each cell within the :cpp:`Box bx`.  If
``USE_CUDA=TRUE``, this incrementation is performed on the GPU, with
GPU optimized loops.

This 4D launch can also be used to work over any sequential set of components, by passing the
number of consecutive components and adding the iterator to the starting component:
:cpp:`fab(i,j,k,n_start+n)`.

The 3D variation of the loop launch does not include a component loop and has the syntax
shown here:

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& fab = mf.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            fab(i,j,k) += 1.;
        });

        /* MACRO VARIATION:
        /
        /   AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        /   {
        /       fab(i,j,k) += 1.;
        /   });
        */
    }

Finally, a 1D version is available for looping over a number of elements, such as particles.
An example of a 1D function launch is given here:

.. highlight:: c++

::

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        Real* AMREX_RESTRICT p = fab.dataPtr();
        const long nitems = fab.box().numPts() * fab.nComp();

        amrex::ParallelFor(nitems,
        [=] AMREX_GPU_DEVICE (long idx)
        {
            p[idx] += 1.;
        });

        /* MACRO VARIATION:
        /
        /   AMREX_PARALLEL_FOR_1D ( nitems, idx,
        /   {
        /       p[idx] += 1.;
        /   });
        */
    }

Instead of passing an :cpp:`Array4`, :cpp:`FArrayBox::dataPtr()` is called to obtain a
CUDA managed pointer to the :cpp:`FArrayBox` data.  This is an alternative way to access
the :cpp:`FArrayBox` data on the GPU. Instead of passing a :cpp:`Box` to define the loop
bounds, a :cpp:`long` or :cpp:`int` number of elements is passed to bound the single
:cpp:`for` loop.  This construct can be used to work on any contiguous set of memory by
passing the number of elements to work on and indexing the pointer to the starting
element: :cpp:`p[idx + 15]`.


Launching general kernels
-------------------------

To launch more general work on the GPU, AMReX provides a standard launch function:
:cpp:`amrex::launch`.  Instead of creating nested loops, this function
prepares the device launch based on a :cpp:`Box`, launches with an appropriate sized
GPU kernel and constructs a thread :cpp:`Box` that defines the work for each thread.
On the CPU, the thread :cpp:`Box` is set equal to the total launch :cpp:`Box`, so
tiling works as expected.  On the GPU, the thread :cpp:`Box` usually
contains a single cell to allow all GPU threads to be utilized effectively.

An example of a generic function launch is shown here:

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);

        amrex::launch(bx,
        [=] AMREX_GPU_DEVICE (Box const& tbx)
        {
            pluseone_array4(tbx, arr);
            FArrayBox fab(arr, tbx.ixType());
            plusone_fab(tbx, fab); // this version takes FArrayBox
        });

        /* MACRO VARIATION
        /
        /   AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        /   {
        /       plusone_array4(tbx, arr);
        /       plusone_fab(tbx, FArrayBox(arr,tbx.ixType()));
        /   });
        */
    }

It also shows how to make a :cpp:`FArrayBox` from :cpp:`Array4` when
needed.  Note that :cpp:`FarrayBox`\ es cannot be passed to GPU
kernels directly.  :cpp:`TilingIfNotGPU()` returns ``false`` in the
GPU case to turn off tiling and maximize the amount of work given to
the GPU in each launch, which substantially improves performance.
When tiling is off, :cpp:`tilebox()` returns the :cpp:`validbox()` of
the :cpp:`FArrayBox` for that iteration.

Offloading work using OpenACC or OpenMP pragmas
-----------------------------------------------

When using OpenACC or OpenMP with AMReX, the GPU offloading work is done
with pragmas placed on the nested loops. This leaves the :cpp:`MFIter` loop
largely unchanged.  An example GPU pragma based :cpp:`MFIter` loop that calls
a Fortran function is given here:

.. highlight:: c++

::

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        plusone_acc(BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_ANYD(fab));
    }

The function ``plusone_acc`` is a CPU host function.  The
:cpp:`FArrayBox` reference
from :cpp:`operator[]` is a reference to a :cpp:`FArrayBox` in host
memory with data that has been placed in managed CUDA memory.
``BL_TO_FORTRAN_BOX`` and ``BL_TO_FORTRAN_ANYD`` behave identically
to implementations used on the CPU.  These macros return the
individual components of the AMReX C++ objects to allow passing to
the Fortran function.

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
unified memory, OpenACC can be told the data is available on the
device using the ``deviceptr`` construct.  For further details
about OpenACC programming, consult the OpenACC user's guide.

The OpenMP implementation of this loop is similar, only requiring
changing the pragmas utilized to obtain the proper offloading. The
OpenMP labelled version of this loop is:

.. highlight:: fortran

::

    !dat = pointer to fab's managed data

    !$omp target teams distribute parallel do collapse(3) schedule(static,1) is_device_ptr(dat)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat(i,j,k) = dat(i,j,k) + 1.0_amrex_real
          end do
       end do
    end do

In this case, ``is_device_ptr`` is used to indicate that :cpp:`dat`
is available in device memory. For further details about programming
with OpenMP for GPU offloading, consult the OpenMP user's guide.


Kernel launch details
---------------------

CUDA kernel calls are asynchronous and they return before the kernel
is finished on the GPU. So the :cpp:`MFIter` loop finishes iterating on
the CPU and is ready to move on to the next work before the actual
work completes on the GPU.  To guarantee consistency,
there is an implicit device synchronization (a GPU barrier) in
the destructor of :cpp:`MFIter`.  This ensures that all GPU work
inside of an :cpp:`MFIter` loop will complete before code outside of
the loop is executed. Any CUDA kernel launches made outside of an
:cpp:`MFIter` loop must ensure appropriate device synchronization
occurs. This can be done by calling :cpp:`Gpu::synchronize()`.

CUDA supports multiple streams and kernels. Kernels launched in the
same stream are executed sequentially, but different streams of kernel
launches may be run in parallel.  For each iteration of :cpp:`MFIter`,
AMReX uses a different CUDA stream (up to 16 streams in total).  This
allows each iteration of an :cpp:`MFIter` loop to run independently,
but in the expected sequence, and maximize the use of GPU parallelism.
However, AMReX uses the default CUDA stream outside of :cpp:`MFIter`
loops.

Launching kernels with AMReX's launch macros or functions implement
a C++ lambda function. Lambdas functions used with CUDA have some
restrictions the user must understand.  First, the function enclosing the
extended lambda must not have private or protected access within its parent
class,  otherwise the code will not compile.  This can be fixed by changing
the access of the enclosing function to public.

Another pitfall that must be considered: if the lambda function
accesses a member of the enclosing class, the lambda function actually
captures :cpp:`this` pointer by value and accesses variables and functions
via :cpp:`this->`.  If the object is not accessible on GPU, the code will
not work as intended.  For example,

.. highlight:: c++

::

    class MyClass {
    public:
        Box bx;
        int m;                           // Unmanaged integer created on the host.
        void f () {
            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                printf("m = %d\n", m);   // Failed attempt to use m on the GPU.
            });
        }
    };

The function ``f`` in the code above will not work unless the :cpp:`MyClass`
object is in unified memory.  If it is undesirable to put the object into
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
            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                printf("m = %d\n", local_m);  // Lambda captures local_m by value.
            });
        }
    };

C++ macros have some important limitations. For example, commas outside
of a set of parentheses are interpreted by the macro, leading to errors such
as:

.. highlight:: c++

::

    AMREX_PARALLEL_FOR_3D (bx, tbx,
    {
        Real a, b;   <---- Error. Macro reads "{ Real a" as a parameter
                                                 and "b; }" as
                                                 another.
        Real a;      <---- OK
        Real b;
    });

One should also avoid using :cpp:`continue` and :cpp:`return` inside the macros
because it is not an actual :cpp:`for` loop.
Users that choose to implement the macro launches should be aware of the limitations
of C++ preprocessing macros to ensure GPU offloading is done properly.

Finally, AMReX's most common CPU threading strategy for GPU/CPU systems is to utilize
OpenMP threads to maintain multi-threaded parallelism on work chosen to run on the host.
This means OpenMP pragmas should be maintained where CPU work is performed and usually
turned off where work is offloaded onto the GPU.  OpenMP pragmas can be turned
off using the conditional pragma and :cpp:`Gpu::notInLaunchRegion()`, as shown below:

.. highlight:: c++

::

    #ifdef AMREX_USE_OMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif

It is generally expected that simply using OpenMP threads to launch GPU work quicker
will show little improvement or even perform worse. So, this conditional statement
should be added to MFIter loops that contain GPU work, unless users specifically test
the performance or are designing more complex workflows that require OpenMP.

.. _sec:gpu:example:

An Example of Migrating to GPU
==============================

The nature of GPU programming poses difficulties for a number
of common AMReX patterns, such as the one below:

.. highlight:: c++

::

   // Given MultiFab uin and uout
   #ifdef AMREX_USE_OMP
   #pragma omp parallel
   #endif
   {
     FArrayBox q;
     for (MFIter mfi(uin,true); mfi.isValid(); ++mfi)
     {
       const Box& tbx = mfi.tilebox();
       const Box& gbx = amrex::grow(tbx,1);
       q.resize(gbx);

       // Do some work with uin[mfi] as input and q as output.
       // The output region is gbx;
       f1(gbx, q, uin[mfi]);

       // Then do more work with q as input and uout[mfi] as output.
       // The output region is tbx.
       f2(tbx, uout[mfi], q);
     }
   }

There are several issues in migrating this code to GPUs that need to
be addressed.  First, functions ``f1`` and ``f2`` have different
work regions (``tbx`` and ``gbx``, respectively) and there are data
dependencies between the two (``q``). This makes it difficult to put
them into a single GPU kernel, so two separate kernels will be
launched, one for each function.

As we have discussed, AMReX uses multiple CUDA streams for launching
kernels.  Because ``q`` is used inside :cpp:`MFIter` loops, multiple
GPU kernels on different streams are accessing its data.  This creates
a race condition.  One way to fix this is to move ``FArrayBox q``
inside the loop to make it local to each loop and use :cpp:`Elixir` to
make it async-safe (see Section :ref:`sec:gpu:classes:elixir`).  This
strategy works well for GPU.  However it is not optimal for OpenMP CPU
threads when CUDA is not used, because of the memory allocation inside
OpenMP parallel region.  It turns out it is actually unnecessary to
make ``FArrayBox q`` local to each iteration when :cpp:`Elixir` is
used to extend the life of its floating point data.  The code below
shows an example of how to rewrite the example in a performance
portable way.

.. highlight:: c++

::

   // Given MultiFab uin and uout
   #ifdef AMREX_USE_OMP
   #pragma omp parallel if (Gpu::notInLaunchRegion())
   #endif
   {
     FArrayBox q;
     for (MFIter mfi(uin,TilingIfNotGPU()); mfi.isValid(); ++mfi)
     {
       const Box& tbx = mfi.tilebox();
       const Box& gbx = amrex::grow(tbx,1);
       q.resize(gbx);
       Elixir eli = q.elixir();
       Array4<Real> const& qarr = q.array();

       Array4<Real const> const& uinarr = uin.const_array(mfi);
       Array4<Real> const& uoutarr = uout.array(mfi);

       amrex::launch(gbx,
       [=] AMREX_GPU_DEVICE (Box const& b)
       {
         f1(b, qarr, uinarr);
       });

       amrex::launch(tbx,
       [=] AMREX_GPU_DEVICE (Box const& b)
       {
         f2(b, uoutarr, qarr);
       });
     }
   }

.. ===================================================================

.. _sec:gpu:assertion:


Assertions, Error Checking and Synchronization
================================================

To help debugging, we often use :cpp:`amrex::Assert` and
:cpp:`amrex::Abort`.  These functions are GPU safe and can be used in
GPU kernels.  However, implementing these functions requires additional
GPU registers, which will reduce overall performance.  Therefore, it
is preferred to implement such calls in debug mode only by wrapping the
calls using ``#ifdef AMREX_DEBUG``.

In CPU code, :cpp:`AMREX_GPU_ERROR_CHECK()` can be called
to check the health of previous GPU launches.  This call
looks up the return message from the most recently completed GPU
launch and aborts if it was not successful. Many kernel
launch macros as well as the :cpp:`MFIter` destructor include a call
to :cpp:`AMREX_GPU_ERROR_CHECK()`. This prevents additional launches
from being called if a previous launch caused an error and ensures
all GPU launches within an :cpp:`MFIter` loop completed successfully
before continuing work.

However, due to asynchronicity, determining the source of the error
can be difficult.  Even if GPU kernels launched earlier in the code
result in a CUDA error, the error may not be output at a nearby call to
:cpp:`AMREX_GPU_ERROR_CHECK()` by the CPU.  When tracking down a CUDA
launch error, :cpp:`Gpu::synchronize()` and
:cpp:`Gpu::streamSynchronize()` can be used to synchronize
the device or the CUDA stream, respectively, and track down the specific
launch that causes the error.

.. ===================================================================


Particle Support
================

.. _sec:gpu:particle:

As with ``MultiFab``, particle data stored in AMReX ``ParticleContainer`` classes are
stored in unified memory when AMReX is compiled with ``USE_CUDA=TRUE``. This means that the :cpp:`dataPtr` associated with particles
is managed and can be passed into GPU kernels. These kernels can be launched with a variety of approaches,
including Cuda C / Fortran and OpenACC. An example Fortran particle subroutine offloaded via OpenACC might
look like the following:

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
       structs(ip)%pos(2) = structs(ip)%pos(2) + uyp(ip)*gaminv(ip)*dt
       structs(ip)%pos(3) = structs(ip)%pos(3) + uzp(ip)*gaminv(ip)*dt
   end do
   !$acc end loop
   !$acc end parallel

   end subroutine push_position_boris

Note the use of the :fortran:`!$acc parallel deviceptr` clause to specify which data has been placed
in managed memory. This instructs OpenACC to treat those variables as if they already live on
the device, bypassing the usual copies. For complete examples of a particle code that has been ported
to GPUs using Cuda, OpenACC, and OpenMP, please see the tutorial `Electromagnetic PIC`_.

.. _`Electromagnetic PIC`: https://amrex-codes.github.io/amrex/tutorials_html/Particles_Tutorial.html#electromagneticpic

GPU-aware implementations of many common particle operations are provided with AMReX, including neighbor list
construction and traversal, particle-mesh deposition and interpolation, parallel reductions of particle data,
and a set of transformation and filtering operations that are useful when operating on sets of particles. For
examples of these features in use, please see :cpp:`Tests/Particles/`.

Finally, the parallel communication of particle data has been ported and optimized for performance on GPU
platforms. This includes :cpp:`Redistribute()`, which moves particles back to the proper grids after their positions
have changed, as well as :cpp:`fillNeighbors()` and :cpp:`updateNeighbors()`, which are used to exchange halo particles.
As with :cpp:`MultiFab` data, these have been designed to minimize host / device traffic as much as possible, and can
take advantage of the Cuda-aware MPI implementations available on platforms such as ORNL's Summit.


Profiling with GPUs
===================

.. _sec:gpu:profiling:

When profiling for GPUs, AMReX recommends ``nvprof``, NVIDIA's visual
profiler.  ``nvprof`` returns data on how long each kernel launch lasted on
the GPU, the number of threads and registers used, the occupancy of the GPU
and recommendations for improving the code.  For more information on how to
use ``nvprof``, see NVIDIA's User's Guide as well as the help web pages of
your favorite supercomputing facility that uses NVIDIA GPUs.

AMReX's internal profilers currently cannot hook into profiling information
on the GPU and an efficient way to time and retrieve that information is
being explored. In the meantime, AMReX's timers can be used to report some
generic timers that are useful in categorizing an application.

Due to the asynchronous launching of GPU kernels, any AMReX timers inside of
asynchronous regions or inside GPU kernels will not measure useful
information.  However, since the :cpp:`MFIter` synchronizes when being
destroyed, any timer wrapped around an :cpp:`MFIter` loop will yield a
consistent timing of the entire set of GPU launches contained within. For
example:

.. highlight:: cpp

::

    BL_PROFILE_VAR("A_NAME", blp);     // Profiling start
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        // gpu works
    }
    BL_PROFILE_STOP(blp);              // Profiling stop

For now, this is the best way to profile GPU codes using ``TinyProfiler``.
If you require further profiling detail, use ``nvprof``.


Performance Tips
================

.. _sec:gpu:performance:

Here are some helpful performance tips to keep in mind when working with
AMReX for GPUs:

* To obtain the best performance when using CUDA kernel launches, all
  device functions called within the launch region should be inlined.
  Inlined functions use substantially fewer registers, freeing up GPU
  resources to perform other tasks. This increases parallel
  performance and greatly reduces runtime.  Functions are written
  inline by putting their definitions in the ``.H`` file and using
  the ``AMREX_FORCE_INLINE`` AMReX macro.  Examples can be found in
  in the `Launch`_ tutorial. For example:

.. highlight:: cpp

::

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void plusone_cudacpp (amrex::Box const& bx, amrex::FArrayBox& fab)
    {
        ...
    }

* Pay attention to what GPUs your job scheduler is assigning to each MPI
  rank. In most cases you'll achieve the best performance when a single
  MPI rank is assigned to each GPU, and has boxes large enough to saturate
  that GPU's compute capacity. While there are some cases where multiple
  MPI ranks per GPU can make sense (typically this would be when you have
  some portion of your code that is not GPU accelerated and want to have
  many MPI ranks to make that part faster), this is probably the minority
  of cases. For example, on OLCF Summit you would want to ensure that your
  resource sets contain one MPI rank and GPU each, using `jsrun -n N -a 1 -c 7 -g 1`,
  where `N` is the total number of MPI ranks/GPUs you want to use. (See the OLCF
  [job step viewer](https://jobstepviewer.olcf.ornl.gov/) for more information.)

  Conversely, if you choose to have multiple GPUs visible to each MPI rank,
  AMReX will attempt to do the best job it can assigning MPI ranks to GPUs by
  doing round robin assignment. This may be suboptimal because this assignment
  scheme would not be aware of locality benefits that come from having an MPI
  rank be on the same socket as the GPU it is managing. If you know the hardware
  layout of the system you're running on, specifically the number of GPUs per
  socket (`M`) and number of GPUs per node (`N`), you can set the preprocessor
  defines `-DAMREX_GPUS_PER_SOCKET=M` and `-DAMREX_GPUS_PER_NODE=N`, which are
  exposed in the GNU Make system through the variables `GPUS_PER_SOCKET` and
  `GPUS_PER_NODE` respectively (see an example in `Tools/GNUMake/sites/Make.olcf`).
  Then AMReX can ensure that each MPI rank selects a GPU on the same socket as
  that rank (assuming your MPI implementation supports MPI 3.)


.. ===================================================================

Inputs Parameters
=================

.. _sec:gpu:parameters:

The following inputs parameters control the behavior of amrex when running on GPUs. They should be prefaced
by "amrex" in your :cpp:`inputs` file.

+----------------------------+-----------------------------------------------------------------------+-------------+-------------+
|                            | Description                                                           |   Type      | Default     |
+============================+=======================================================================+=============+=============+
| use_gpu_aware_mpi          | Whether to use GPU memory for communication buffers during MPI calls. | Bool        | False       |
|                            | If true, the buffers will use device memory. If false, they will use  |             |             |
|                            | pinned memory. In practice, we find it is usually not worth it to use |             |             |
|                            | GPU aware MPI.                                                        |             |             |
+----------------------------+-----------------------------------------------------------------------+-------------+-------------+
| abort_on_out_of_gpu_memory | If the size of free memory on the GPU is less than the size of a      | Bool        | False       |
|                            | requested allocation, AMReX will call AMReX::Abort() with an error    |             |             |
|                            | describing how much free memory there is and what was requested.      |             |             |
+----------------------------+-----------------------------------------------------------------------+-------------+-------------+

Basic Gpu Debugging
===================

- Turn off GPU offloading for some part of the code with

.. highlight:: cpp

::

    Gpu::setLaunchRegion(0);
    ... ;
    Gpu::setLaunchRegion(1);

Note that functions, ``amrex::launch`` and ``amrex::ParallelFor``, do
not respect the launch region flag.  Only the macros (e.g.,
``AMREX_LAUNCH_HOST_DEVICE_LAMBDA`` and ``AMREX_HOST_DEVICE_FOR_*D``) do.

Cuda-specific tests
-------------------

- To test if your kernels have launched, run

::

    nvprof ./main3d.xxx

- Run under ``nvprof -o profile%p.nvvp ./main3d.xxxx`` for
  a small problem and examine page faults using nvvp

- Run under ``cuda-memcheck``

- Run under ``cuda-gdb``

- Run with ``CUDA_LAUNCH_BLOCKING=1``.  This means that only one
  kernel will run at a time.  This can help identify if there are race
  conditions.
