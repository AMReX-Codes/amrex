include_guard(GLOBAL)

#
# Resolve the target CUDA architecture(s) into the standard CMAKE_CUDA_ARCHITECTURES
# cache variable.
#
# This must be included BEFORE enable_language(CUDA) so that the CMP0104 NEW behavior
# does not initialize CMAKE_CUDA_ARCHITECTURES to a compiler-default architecture, which
# would defeat both our "native" default and the precedence of the legacy AMReX hints.
#
# We standardize on the CMake-native interface:
#   - CMAKE_CUDA_ARCHITECTURES cache variable (e.g. -DCMAKE_CUDA_ARCHITECTURES=80)
#   - CUDAARCHS environment variable (CMake's initializer for the above)
#
# The legacy AMReX hints are still honored (with a deprecation warning), keeping their
# historical precedence:
#   0. -DCMAKE_CUDA_ARCHITECTURES=... given on the command line of this configure step
#   1. AMReX_CUDA_ARCH cache variable  (deprecated -> CMAKE_CUDA_ARCHITECTURES)
#   2. AMREX_CUDA_ARCH  environment variable (deprecated -> CUDAARCHS)
#   3. CMAKE_CUDA_ARCHITECTURES cache variable
#   4. CUDAARCHS environment variable
#   5. "native" (build for the locally installed GPU)
#

set(_amrex_cuda_archs_legacy FALSE)

#
# Was CMAKE_CUDA_ARCHITECTURES chosen by the user (or by the parent project), or is it
# merely the value CMake initialized from the compiler default? The latter happens when
# AMReX is used as a subproject and the parent enables the CUDA language before adding
# AMReX: CMP0104 NEW then already stored nvcc's default architecture in the cache (e.g.
# sm_52 with CUDA 12.x, below the compute capability AMReX supports), which would silently
# win over AMReX's own default. CMake writes that entry with its own help string, so the
# two cases can be told apart; a value that came from CUDAARCHS is picked up below anyway.
#
set(_amrex_cuda_archs_given FALSE)
set(_amrex_cuda_archs_cmdline FALSE)
if (DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(_amrex_cuda_archs_given TRUE)
   get_property(_amrex_cuda_archs_help CACHE CMAKE_CUDA_ARCHITECTURES PROPERTY HELPSTRING)
   if (_amrex_cuda_archs_help MATCHES "specified on the command line")
      # -DCMAKE_CUDA_ARCHITECTURES=... for this configure step: more explicit than the
      # deprecated AMReX hints, which may just be stale entries of an older build directory
      set(_amrex_cuda_archs_cmdline TRUE)
   elseif (CMAKE_CUDA_COMPILER_LOADED AND _amrex_cuda_archs_help STREQUAL "CUDA architectures")
      # CMake's own help string: either the compiler default stored by enable_language(CUDA)
      # - which happens when a parent project enables CUDA before adding AMReX, and which
      # AMReX must not mistake for a deliberate choice - or the value of CUDAARCHS, which
      # is honored below anyway. A parent that set the variable after enable_language(CUDA)
      # is recognized by its value differing from the cache entry.
      get_property(_amrex_cuda_archs_cached CACHE CMAKE_CUDA_ARCHITECTURES PROPERTY VALUE)
      if (CMAKE_CUDA_ARCHITECTURES STREQUAL _amrex_cuda_archs_cached
          AND NOT CMAKE_CUDA_ARCHITECTURES STREQUAL "$ENV{CUDAARCHS}")
         set(_amrex_cuda_archs_given FALSE)
         message(STATUS
            "The CUDA language was enabled before AMReX was added, so CMake preset "
            "CMAKE_CUDA_ARCHITECTURES to the compiler default (${CMAKE_CUDA_ARCHITECTURES}); "
            "applying AMReX's own architecture selection instead. Set "
            "CMAKE_CUDA_ARCHITECTURES (or CUDAARCHS) before enable_language(CUDA) to "
            "choose the architectures yourself.")
      endif ()
      unset(_amrex_cuda_archs_cached)
   endif ()
   unset(_amrex_cuda_archs_help)
endif ()

if (_amrex_cuda_archs_cmdline AND (DEFINED AMReX_CUDA_ARCH OR DEFINED ENV{AMREX_CUDA_ARCH}))
   message(WARNING
      "Both CMAKE_CUDA_ARCHITECTURES and the deprecated AMReX_CUDA_ARCH/AMREX_CUDA_ARCH "
      "were given; using the CMAKE_CUDA_ARCHITECTURES value passed on the command line "
      "(${CMAKE_CUDA_ARCHITECTURES}).")
   set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES}")
   # a leftover entry of a build directory configured with an older AMReX must not come
   # back on the next configure step
   unset(AMReX_CUDA_ARCH CACHE)
elseif (DEFINED AMReX_CUDA_ARCH)
   message(WARNING
      "AMReX_CUDA_ARCH is deprecated for CMake builds; set the standard "
      "CMAKE_CUDA_ARCHITECTURES instead (e.g. -DCMAKE_CUDA_ARCHITECTURES=80, or 'native'). "
      "Its value takes precedence over CMAKE_CUDA_ARCHITECTURES for this configure step "
      "and is then removed from the cache.")
   set(_amrex_cuda_archs "${AMReX_CUDA_ARCH}")
   set(_amrex_cuda_archs_legacy TRUE)
   # Do not leave a stale cache entry behind: it would silently keep overriding
   # -DCMAKE_CUDA_ARCHITECTURES on every subsequent configure of this build directory,
   # without an option of that name being visible in the cache anymore.
   unset(AMReX_CUDA_ARCH CACHE)
elseif (DEFINED ENV{AMREX_CUDA_ARCH})
   message(WARNING
      "The AMREX_CUDA_ARCH environment variable is deprecated for CMake builds (it remains "
      "in use for GNU Make builds); use the standard CUDAARCHS environment variable instead "
      "(or -DCMAKE_CUDA_ARCHITECTURES=...).")
   set(_amrex_cuda_archs "$ENV{AMREX_CUDA_ARCH}")
   set(_amrex_cuda_archs_legacy TRUE)
elseif (_amrex_cuda_archs_given)
   set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES}")
elseif (DEFINED ENV{CUDAARCHS})
   set(_amrex_cuda_archs "$ENV{CUDAARCHS}")
else ()
   set(_amrex_cuda_archs "native")
endif ()

#
# Back-compatibility normalization, applied to the legacy hints only:
#   - "Auto" (the historical default) -> "native"
#   - "All" -> "all" and "Common" -> "all-major"
#   - legacy NVIDIA generation names to their compute capability (e.g. Volta -> 70), as the
#     deprecated FindCUDA helper used to accept via AMReX_CUDA_ARCH, with a "+Tegra" or
#     "+Tesla" suffix dropped. GPU SASS code is forward-compatible
#     across minor revisions of the same generation, so the major base value covers
#     the whole family (e.g. 80 runs on 86/87). "Blackwell" is the exception: it
#     spans two binary-incompatible families (data-center sm_100 and consumer
#     sm_120), so it expands to both. The generations below AMReX's minimum are
#     translated as well, only so that the compute capability check below can report them
#     by number instead of letting CMake fail on an unrecognized string.
#   - drop the "+PTX" suffix ("7.0+PTX" -> "70"): the plain integer form of
#     CUDA_ARCHITECTURES already embeds PTX in addition to the SASS code
#   - strip the decimal dot per entry ("8.0" -> "80", "7.5" -> "75")
# Values given through the modern interfaces are passed through unchanged, i.e. everything
# the CUDA_ARCHITECTURES target property understands:
#   native / all / all-major / integers / <NN>a / -real / -virtual suffixes.
#
if (_amrex_cuda_archs_legacy)
   set(_amrex_cuda_archs_norm)
   foreach (_arch IN LISTS _amrex_cuda_archs)
      # all values CUDA_ARCHITECTURES accepts are lower case
      string(TOLOWER "${_arch}" _arch)
      string(REGEX REPLACE "\\+ptx$" "" _arch "${_arch}")
      string(REGEX REPLACE "\\+(tegra|tesla)$" "" _arch "${_arch}")
      if (_arch STREQUAL "auto")
         set(_arch "native")
      elseif (_arch STREQUAL "common")
         set(_arch "all-major")
      elseif (_arch STREQUAL "fermi")
         set(_arch "20")
      elseif (_arch STREQUAL "kepler")
         set(_arch "35")
      elseif (_arch STREQUAL "maxwell")
         set(_arch "50")
      elseif (_arch STREQUAL "pascal")
         set(_arch "60")
      elseif (_arch STREQUAL "volta")
         set(_arch "70")
      elseif (_arch STREQUAL "turing")
         set(_arch "75")
      elseif (_arch STREQUAL "ampere")
         set(_arch "80")
      elseif (_arch STREQUAL "ada" OR _arch STREQUAL "lovelace")
         set(_arch "89")
      elseif (_arch STREQUAL "hopper")
         set(_arch "90")
      elseif (_arch STREQUAL "blackwell")
         set(_arch "100;120")
      else ()
         string(REPLACE "." "" _arch "${_arch}")
      endif ()
      list(APPEND _amrex_cuda_archs_norm "${_arch}")
   endforeach ()
   unset(_arch)
   set(_amrex_cuda_archs "${_amrex_cuda_archs_norm}")
   unset(_amrex_cuda_archs_norm)
endif ()

# the aliases native/all/all-major are only resolved once the CUDA language is enabled,
# so those are checked in AMReXCUDAOptions instead
amrex_filter_cuda_archs(_amrex_cuda_archs FALSE)

set(CMAKE_CUDA_ARCHITECTURES "${_amrex_cuda_archs}" CACHE STRING
   "CUDA architectures: 'native', 'all-major', or e.g. 80;90a (see CMake CUDA_ARCHITECTURES)" FORCE)

unset(_amrex_cuda_archs)
unset(_amrex_cuda_archs_legacy)
unset(_amrex_cuda_archs_given)
unset(_amrex_cuda_archs_cmdline)
