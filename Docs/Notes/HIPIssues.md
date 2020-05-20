
# Critical

* HIP HCC has a will-not-fix bug in aggregate initialization.  This
  results in compilation failure in many cases.  This is supposed to
  be fixed in HIP Clang.

* Math functions are currently not in std namespace.  We have added
  some math functions to `namespace amrex::Math`, but we do not want
  to add all functions in `cmath` there.  Sometimes using `std` math
  functions results in compilation failure, whereas sometimes it
  results in runtime error with a cryptic error message.  This is
  supposed to be fixed in HIP Clang.

# Major

* `printf` does not work in device functions.

