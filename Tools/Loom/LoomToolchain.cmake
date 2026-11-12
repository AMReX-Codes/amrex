# CMake toolchain for AMReX_GPU_BACKEND=LOOM.
#
# C++ sources compile with Loom's single-source offload driver, loom-cxx: it
# parses each source once with cxx and emits one object holding the host code,
# the device code of the launched kernels and the kernel descriptors that host
# code references. Executables and shared libraries link with clang++ against
# libc++, the C++ library the host code was compiled for. AMReX adds the device
# options (target, code object registration, drop policy) and libhrx; see
# Tools/CMake/AMReXParallelBackends.cmake.
#
# Configure with -DCMAKE_TOOLCHAIN_FILE=<amrex>/Tools/Loom/LoomToolchain.cmake
# and these settings, as -D cache variables or environment variables (see
# Tools/Loom/single-source/env.sh):
#   AMREX_LOOM_CXX_DRIVER       loom-cxx
#   AMREX_LOOM_LINK_CXX         clang++ that links
#   AMREX_LOOM_LIBCXX_ROOT      libc++ install (include/c++/v1 and lib)
#   AMREX_LOOM_GCC_INSTALL_DIR  GCC's lib/gcc/<triple>/<version>: compiler
#                               headers for compiles, crt files for links
# Optional:
#   AMREX_LOOM_GLIBC_INCLUDE    glibc headers (default /usr/include)
#   AMREX_LOOM_SYSROOT          sysroot holding glibc's libraries for links
#   AMREX_LOOM_HOST_TRIPLE      host triple (default x86_64-linux-gnu)

foreach(_var AMREX_LOOM_CXX_DRIVER AMREX_LOOM_LINK_CXX AMREX_LOOM_LIBCXX_ROOT
             AMREX_LOOM_GCC_INSTALL_DIR AMREX_LOOM_GLIBC_INCLUDE
             AMREX_LOOM_SYSROOT AMREX_LOOM_HOST_TRIPLE)
  if(NOT DEFINED ${_var} AND DEFINED ENV{${_var}})
    set(${_var} "$ENV{${_var}}")
  endif()
  # try_compile projects (compiler checks) see the same settings.
  list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES ${_var})
endforeach()
foreach(_var AMREX_LOOM_CXX_DRIVER AMREX_LOOM_LINK_CXX AMREX_LOOM_LIBCXX_ROOT
             AMREX_LOOM_GCC_INSTALL_DIR)
  if(NOT ${_var})
    message(FATAL_ERROR "LoomToolchain.cmake: set ${_var}")
  endif()
endforeach()
if(NOT AMREX_LOOM_GLIBC_INCLUDE)
  set(AMREX_LOOM_GLIBC_INCLUDE /usr/include)
endif()
if(NOT AMREX_LOOM_HOST_TRIPLE)
  set(AMREX_LOOM_HOST_TRIPLE x86_64-linux-gnu)
endif()

set(CMAKE_CXX_COMPILER "${AMREX_LOOM_CXX_DRIVER}")
# loom-cxx discovers no include directories: libc++, then the compiler's
# headers, then libc, as a clang++ -stdlib=libc++ compile searches them.
# Without -std, cxx parses C++26, and libc++'s C++26 headers do not parse yet;
# start from C++20, which AMReX requires. CMake identifies this default and
# adds a later -std for targets that need a newer standard.
string(JOIN " " CMAKE_CXX_FLAGS_INIT
       "--target=${AMREX_LOOM_HOST_TRIPLE} -std=c++20"
       "-isystem ${AMREX_LOOM_LIBCXX_ROOT}/include/c++/v1"
       "-isystem ${AMREX_LOOM_GCC_INSTALL_DIR}/include"
       "-isystem ${AMREX_LOOM_GLIBC_INCLUDE}")

# Link with clang++ and libc++. The compile flags above are loom-cxx's and
# stay off the link line. cxx emits non-PIC code, so executables link without
# PIE.
set(_amrex_loom_link "\"${AMREX_LOOM_LINK_CXX}\"")
string(APPEND _amrex_loom_link
       " --gcc-install-dir=${AMREX_LOOM_GCC_INSTALL_DIR}")
if(AMREX_LOOM_SYSROOT)
  string(APPEND _amrex_loom_link " --sysroot=${AMREX_LOOM_SYSROOT}")
endif()
string(APPEND _amrex_loom_link
       " -stdlib=libc++ -L${AMREX_LOOM_LIBCXX_ROOT}/lib"
       " -Wl,-rpath,${AMREX_LOOM_LIBCXX_ROOT}/lib")
string(JOIN " " CMAKE_CXX_LINK_EXECUTABLE
       "${_amrex_loom_link} -no-pie <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS>"
       "<OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
string(JOIN " " CMAKE_CXX_CREATE_SHARED_LIBRARY
       "${_amrex_loom_link} <CMAKE_SHARED_LIBRARY_CXX_FLAGS> <LINK_FLAGS>"
       "<CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> <SONAME_FLAG><TARGET_SONAME>"
       "-o <TARGET> <OBJECTS> <LINK_LIBRARIES>")
set(CMAKE_CXX_CREATE_SHARED_MODULE "${CMAKE_CXX_CREATE_SHARED_LIBRARY}")
unset(_amrex_loom_link)
