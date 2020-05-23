
# Critical

* Global variables.  Could DPC++ support global variables and add
  something similar to cudaMemcpyToSymbol?
  [oneAPI-spec issue #125](https://github.com/oneapi-src/oneAPI-spec/issues/125)

* Device API for random number generator.  Currently we can only use
  oneMKL's host API to generate random numbers.
  [Intel oneAPI Base Toolkit Forum](https://software.intel.com/en-us/forums/intel-oneapi-base-toolkit/topic/856436)
  [oneAPI-spec issue #139](https://github.com/oneapi-src/oneAPI-spec/issues/139)

* Recursive function call on device.  This is very important for ECP
  WarpX code.
  [oneAPI-spec issue #123](https://github.com/oneapi-src/oneAPI-spec/issues/123)

* Memory fence.  Could DPC++ privode a memory fence function for the
  whole device (not just group)?  Or is the CUDA distinction between
  `__threadfence` and `__thread_block` unnecessary for Intel GPUs?
  [oneAPI-spec issue #130](https://github.com/oneapi-src/oneAPI-spec/issues/130)

* The compiler has some troubles with some very big device functions
  (e.g., `mlndlap_stencil_rap` in
  `Src/LinearSolvers/MLMG/AMReX_MLNodeLap_3D_K.H`).  It hangs at JIT
  compilation.  We have to disable GPU launch for these functions and
  run them on CPU. 

# Major

* The maximum size of kernel parameter is 1KB on current Intel GPUs.
  This is not sufficient for many of our kernels.

* Sometimes the JIT compilation will raise floating-point exception at
  runtime.  This forces us to disable floating-point exception signal
  handling that we often rely on for debugging.

* Option to be less OOP.  Could we have access to thread id, group id,
  memory fence, barrier functions, etc. without using an nd_item like
  object?
  [oneAPI-spec issue #118](https://github.com/oneapi-src/oneAPI-spec/issues/118)

* Local memory.  Could DPC++ support static local memory
  (e.g. something like CUDA `__shared__ a[256]`) and dynamic local
  memory (e.g., something like CUDA `extern __shared__ a[]` with the
  amount of memory specified at runtime during kernel launch) from
  anywhere in device code?
  [oneAPI-spec issue #126](https://github.com/oneapi-src/oneAPI-spec/issues/126)

# Minor

* Compiler flag to make implicit capture of this pointer via `[=]` an
  error.  [Implicit capture of this pointer](http://eel.is/c++draft/depr#capture.this)
  has been deprecated in C++ 20.  For many codes, it's almost always a
  bug when `this` is implicitly captured onto device via `[=]`.
  [oneAPI-spec issue #127](https://github.com/oneapi-src/oneAPI-spec/issues/127)

* Host callback.  Could DPC++ support appending a host callback
  function to an ordered queue?
  [oneAPI-spec issue #124](https://github.com/oneapi-src/oneAPI-spec/issues/124)

* Subgroup size.  Querying `sycl::info::device::sub_group_size` gives
  several numbers.  For example, we get 8, 16 and 32 for Gen9.  We
  would like to specify the sub group size and this feature is
  supported.  All three sizes seem to work except that subgroup
  primitives such as `shuffle_down` do not work for all sizes.  By try
  and error, we have found that shuffle_down works for 16.  Could
  oneAPI provide a query function for returning the "primitive"
  subgroup size?
  [oneAPI-spec issue #118](https://github.com/oneapi-src/oneAPI-spec/issues/118)

* `assert(0)`. `assert(0)` when called on device does not throw any
  errors or abort the run.  Is it possible to make it abort or return
  an error code that can be checked on the host?  In CUDA, the users
  can check an error code.
  [oneAPI-spec issue #128](https://github.com/oneapi-src/oneAPI-spec/issues/128)

* `sycl::abs`. `sycl::abs(int)` returns an `unsigned int` in contrast to
  `int std::abs(int)`.  Currently `std::abs` does not work on device.  If
  `std::abs` is made to work on device, could it return the same type
  as the C++ standard?
  [oneAPI-spec issue #129](https://github.com/oneapi-src/oneAPI-spec/issues/129)

# Resolved

* ~~Classes that are not standard layout.  The current specification of
  oneAPI does not support the capture of objects that are not standard
  layout.  This includes the following example,~~

  ```
  class A {int a;}; class B {long B;}; class C : A, B {};
  ```

  ~~AMReX has a data structure called GpuTuple that is built with a
  pattern like the example shown above.  It works in CUDA, but not in
  DPC++.  We wish this requirement can be relaxed.~~

  This restriction has been relaxed since beta5.
