
* Recursive function call on device.  This is very important for ECP
  WarpX code.

* Subgrourp size.  Querying `sycl::info::device::sub_group_size` gives
  several numbers.  For example, we get 8, 16 and 32 for Gen9.  We
  would like to specify the sub group size and this feature is
  supported.  All three sizes seem to work except that subgroup
  primitives such as `shuffle_down` do not work for all sizes.  By try
  and error, we have found that shuffle_down works for 16.  Could
  oneAPI provide a query function for returning the "primitive"
  subgroup size?

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

* Host callback.  Could DPC++ support appending a host callback
  function to an ordered queue?

* Global variables.  Could DPC++ support global variables and add
  something similar to cudaMemcpyToSymbol?

* Option to be less OOP.  Could we have access to thread id, group id,
  memory fence, barrier functions, etc. without using an nd_item like
  object?

* Local memory.  Could DPC++ support static local memory
  (e.g. something like CUDA `__shared__ a[256]`) and dynamic local
  memory (e.g., something like CUDA `extern __shared__ a[]` with the
  amount of memory specified at runtime during kernel launch) from
  anywhere in device code?

* Compiler flag to make implicit capture of this pointer via `[=]` an
  error.  [Implicit capture of this pointer](http://eel.is/c++draft/depr#capture.this)
  has been deprecated in C++ 20.  For many codes, it's almost always a
  bug when `this` is implicitly captured onto device via `[=]`.

* `assert(0)`. `assert(0)` when called on device does not throw any
  errors or abort the run.  Is it possible to make it abort or return
  an error code that can be checked on the host?  In CUDA, the users
  can check an error code.

* `sycl::abs`. `sycl::abs(int)` returns an `unsigned int` in contrast to
  `int std::abs(int)`.  Currently `std::abs` does not work on device.  If
  `std::abs` is made to work on device, could it return the same type
  as the C++ standard?

* Memory fence.  Could DPC++ privode a memory fence function for the
  whole device (not just group)?  Or is the CUDA distinction between
  `__threadfence` and `__thread_block` unnecessary for Intel GPUs?
