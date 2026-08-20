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
#   0. -DCMAKE_CUDA_ARCHITECTURES=... given on the command line (and remembered for the
#      automatic re-configure steps of that build directory)
#   1. AMReX_CUDA_ARCH cache variable  (deprecated -> CMAKE_CUDA_ARCHITECTURES)
#   2. AMREX_CUDA_ARCH  environment variable (deprecated -> CUDAARCHS)
#   3. CMAKE_CUDA_ARCHITECTURES cache variable
#   4. CUDAARCHS environment variable
#   5. "native" (build for the locally installed GPU)
#

set(_amrex_cuda_archs_legacy FALSE)

#
# Is a normal (non-cache) variable of this name in effect, i.e. one that shadows the cache
# entry? It is honored like any other user choice below, but it also reaches places AMReX
# cannot correct, so remember it for the diagnostic at the end of this file.
#
set(_amrex_cuda_archs_shadow "")
set(_amrex_cuda_archs_shadowed FALSE)
set(_amrex_cuda_archs_cached "")
# whether there is a cache entry at all, which its value cannot answer: CMake also takes the
# false values (CMAKE_CUDA_ARCHITECTURES=OFF, for packagers who pass the flags themselves)
set(_amrex_cuda_archs_had_cache FALSE)
if (DEFINED CMAKE_CUDA_ARCHITECTURES)
   if (NOT DEFINED CACHE{CMAKE_CUDA_ARCHITECTURES})
      set(_amrex_cuda_archs_shadowed TRUE)
   else ()
      set(_amrex_cuda_archs_had_cache TRUE)
      get_property(_amrex_cuda_archs_cached CACHE CMAKE_CUDA_ARCHITECTURES PROPERTY VALUE)
      if (NOT CMAKE_CUDA_ARCHITECTURES STREQUAL _amrex_cuda_archs_cached)
         set(_amrex_cuda_archs_shadowed TRUE)
      endif ()
   endif ()
   if (_amrex_cuda_archs_shadowed)
      set(_amrex_cuda_archs_shadow "${CMAKE_CUDA_ARCHITECTURES}")
   endif ()
endif ()

#
# Hints that are set but blank ("export AMREX_CUDA_ARCH=", -DAMReX_CUDA_ARCH=, or nothing
# but spaces) are treated as not given: they satisfy DEFINED, but an empty architecture list
# is rejected by CMake instead of falling back to the default selection below - and the blank
# would only be recognized as empty after the normalization further down, once it had already
# been chosen over every other hint.
#
if (DEFINED AMReX_CUDA_ARCH AND AMReX_CUDA_ARCH MATCHES "^[ \t\r\n]*$")
   unset(AMReX_CUDA_ARCH CACHE)
   unset(AMReX_CUDA_ARCH)
endif ()
if (DEFINED ENV{AMREX_CUDA_ARCH} AND "$ENV{AMREX_CUDA_ARCH}" MATCHES "^[ \t\r\n]*$")
   unset(ENV{AMREX_CUDA_ARCH})
endif ()
if (DEFINED ENV{CUDAARCHS} AND "$ENV{CUDAARCHS}" MATCHES "^[ \t\r\n]*$")
   unset(ENV{CUDAARCHS})
endif ()

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
set(_amrex_cuda_archs_cmdline_now FALSE)

# -DAMReX_CUDA_ARCH=... passed to this configure step, as opposed to a leftover cache entry
# of an earlier one: it outranks a command line choice this build directory only remembers
set(_amrex_cuda_arch_cmdline_now FALSE)
if (DEFINED CACHE{AMReX_CUDA_ARCH})
   get_property(_amrex_cuda_arch_help CACHE AMReX_CUDA_ARCH PROPERTY HELPSTRING)
   if (_amrex_cuda_arch_help MATCHES "specified on the command line")
      set(_amrex_cuda_arch_cmdline_now TRUE)
   endif ()
   unset(_amrex_cuda_arch_help)
endif ()

if (DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(_amrex_cuda_archs_given TRUE)
   get_property(_amrex_cuda_archs_help CACHE CMAKE_CUDA_ARCHITECTURES PROPERTY HELPSTRING)
   if (_amrex_cuda_archs_help MATCHES "specified on the command line")
      # -DCMAKE_CUDA_ARCHITECTURES=... for this configure step: more explicit than the
      # deprecated AMReX hints, which may just be stale entries of an older build directory
      set(_amrex_cuda_archs_cmdline TRUE)
      set(_amrex_cuda_archs_cmdline_now TRUE)
   elseif (DEFINED AMREX_CUDA_ARCHS_CMDLINE
           AND CMAKE_CUDA_ARCHITECTURES STREQUAL AMREX_CUDA_ARCHS_CMDLINE)
      # CMake marks an entry as "specified on the command line" only for the configure step
      # -D... is passed to, and writing the resolved architectures back to the cache below
      # replaces that help string anyway. Keep honoring the choice on the automatic
      # re-configure steps that follow (a build system regenerates itself when a
      # CMakeLists.txt changes), as long as the cache entry still holds what the command
      # line resolved to. Passing a different -D value updates the marker, changing the
      # entry by other means drops it, both at the end of this file. A deprecated hint that
      # is itself passed to this configure step is newer than what is remembered here and
      # takes over again (see the next block).
      set(_amrex_cuda_archs_cmdline TRUE)
   elseif (CMAKE_CUDA_COMPILER_LOADED AND _amrex_cuda_archs_help STREQUAL "CUDA architectures")
      # CMake's own help string: either the compiler default stored by enable_language(CUDA)
      # - which happens when a parent project enables CUDA before adding AMReX, and which
      # AMReX must not mistake for a deliberate choice - or the value of CUDAARCHS, which
      # is honored below anyway. A parent that set the variable after enable_language(CUDA)
      # is recognized by its value differing from the cache entry.
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
   endif ()
   unset(_amrex_cuda_archs_help)
endif ()

if (_amrex_cuda_archs_cmdline AND (_amrex_cuda_archs_cmdline_now OR NOT _amrex_cuda_arch_cmdline_now)
    AND (DEFINED AMReX_CUDA_ARCH OR DEFINED ENV{AMREX_CUDA_ARCH}))
   # a normal variable shadowing the cache entry supplies the value that is in effect, which
   # is then not the one the command line carries (reported by the diagnostic further down)
   if (_amrex_cuda_archs_shadowed)
      set(_amrex_cuda_archs_source "the CMAKE_CUDA_ARCHITECTURES value in effect")
   elseif (_amrex_cuda_archs_cmdline_now)
      set(_amrex_cuda_archs_source "the CMAKE_CUDA_ARCHITECTURES value from the command line")
   else ()
      set(_amrex_cuda_archs_source
          "the CMAKE_CUDA_ARCHITECTURES value this build directory remembers from its command line")
   endif ()
   message(WARNING
      "Both CMAKE_CUDA_ARCHITECTURES and the deprecated AMReX_CUDA_ARCH/AMREX_CUDA_ARCH "
      "were given; using ${_amrex_cuda_archs_source} (${CMAKE_CUDA_ARCHITECTURES}).")
   unset(_amrex_cuda_archs_source)
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
   # whatever the command line of an earlier configure step chose has just been replaced
   set(_amrex_cuda_archs_cmdline FALSE)
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
   set(_amrex_cuda_archs_cmdline FALSE)
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
   # GNU Make builds take a whitespace-separated list (AMREX_CUDA_ARCH="70 80"), which CMake
   # would see as a single architecture named "70 80"
   string(STRIP "${_amrex_cuda_archs}" _amrex_cuda_archs)
   string(REGEX REPLACE "[ \t\r\n]+" ";" _amrex_cuda_archs "${_amrex_cuda_archs}")

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

# The architectures are selected again on every configure step, so a hint that is still in
# place can replace what an earlier step chose: the deprecated AMREX_CUDA_ARCH environment
# variable outranks the CMAKE_CUDA_ARCHITECTURES cache entry, for instance, and is honored
# again once the -DAMReX_CUDA_ARCH=... of that earlier step is gone. Nothing on this command
# line asked for the change, so do not let the build directory change quietly.
#
# Only a step that follows one of AMReX's own can tell: on a first configure the cache entry
# holds whatever the caller supplied with it - an initial-cache script (-C), a toolchain file
# or a parent project - and resolving that into the architectures AMReX builds for is not a
# change of anything. The architectures a parent project preset before adding AMReX are
# reported by their own message above.
if (DEFINED CACHE{AMREX_CUDA_ARCHS_CONFIGURED} AND _amrex_cuda_archs_had_cache
    AND _amrex_cuda_archs_given
    AND NOT _amrex_cuda_archs STREQUAL _amrex_cuda_archs_cached
    AND NOT _amrex_cuda_archs_cmdline_now AND NOT _amrex_cuda_arch_cmdline_now)
   message(WARNING
      "The CUDA architectures of this build directory change from "
      "'${_amrex_cuda_archs_cached}' to '${_amrex_cuda_archs}', although none were requested "
      "on this command line: they are selected again on every configure step, from the hints "
      "that are still in place (see the messages above). Pass "
      "-DCMAKE_CUDA_ARCHITECTURES=<archs> to pin the choice for this build directory.")
endif ()

# A normal (non-cache) variable shadows the cache entry, so it also hides a
# -DCMAKE_CUDA_ARCHITECTURES=... of this configure step, which never becomes visible.
if (_amrex_cuda_archs_shadowed AND _amrex_cuda_archs_cmdline_now)
   message(WARNING
      "CMAKE_CUDA_ARCHITECTURES was passed on the command line "
      "(${_amrex_cuda_archs_cached}), but a normal (non-cache) variable of that name, set "
      "by a parent project or by a toolchain file, holds '${_amrex_cuda_archs_shadow}' and "
      "shadows it, so AMReX is built for '${_amrex_cuda_archs}'. Drop that assignment or "
      "write it to the cache - set(CMAKE_CUDA_ARCHITECTURES <archs> CACHE STRING \"\") - "
      "for the command line to take effect.")
endif ()

if (_amrex_cuda_archs_shadowed AND NOT _amrex_cuda_archs_shadow STREQUAL _amrex_cuda_archs)
   if (CMAKE_TOOLCHAIN_FILE)
      # A toolchain file is read again in every project CMake configures on the side, in
      # particular in the compiler test projects of enable_language(CUDA) below. An
      # unconditional set(CMAKE_CUDA_ARCHITECTURES ...) in such a file therefore restores
      # its own value there, after the value CMake passes in, and AMReX cannot correct it -
      # unlike the same assignment in a parent project, which the normal variable below
      # takes care of.
      message(WARNING
         "CMAKE_CUDA_ARCHITECTURES is set as a normal (non-cache) variable holding "
         "'${_amrex_cuda_archs_shadow}', while AMReX is built for '${_amrex_cuda_archs}'. "
         "If that assignment comes from the toolchain file (${CMAKE_TOOLCHAIN_FILE}), it "
         "also applies to the compiler test projects CMake configures next, which read the "
         "toolchain file again: those fail with e.g. \"nvcc fatal : Unsupported gpu "
         "architecture\" whenever the CUDA compiler cannot build one of the requested "
         "architectures, even though AMReX itself would use the value above. Request "
         "architectures the compiler supports, or write them as a cache entry - "
         "set(CMAKE_CUDA_ARCHITECTURES <archs> CACHE STRING \"\") - which lets AMReX's "
         "selection take effect everywhere.")
   else ()
      # the enclosing project keeps its own value for the targets it adds itself
      message(STATUS
         "CMAKE_CUDA_ARCHITECTURES is set as a normal (non-cache) variable holding "
         "'${_amrex_cuda_archs_shadow}' in the enclosing scope; AMReX is built for "
         "${_amrex_cuda_archs}, while targets added outside of it keep the former.")
   endif ()
endif ()
unset(_amrex_cuda_archs_shadow)

# A normal variable of the same name - set by a parent project or by a toolchain file -
# survives the cache write above (CMP0126 is NEW for the CMake 3.25 required here) and
# shadows it for the rest of this scope, i.e. for enable_language(CUDA) and every target
# added afterwards: the architectures resolved here would be reported and exported, but the
# build would use the unfiltered value. Let the resolved list win in AMReX's scope.
set(CMAKE_CUDA_ARCHITECTURES "${_amrex_cuda_archs}")

# this build directory has been through the selection above, which the next configure step
# needs to know to tell a changed selection from the value a first configure was handed
set(AMREX_CUDA_ARCHS_CONFIGURED TRUE CACHE INTERNAL
   "AMReX has resolved the CUDA architectures of this build directory")

# see the "specified on the command line" check above; a value that a shadowing normal
# variable supplied did not come from the command line and must not be remembered as if it did
if (_amrex_cuda_archs_cmdline AND NOT _amrex_cuda_archs_shadowed)
   set(AMREX_CUDA_ARCHS_CMDLINE "${_amrex_cuda_archs}" CACHE INTERNAL
      "CMAKE_CUDA_ARCHITECTURES as it was resolved from the command line")
else ()
   unset(AMREX_CUDA_ARCHS_CMDLINE CACHE)
endif ()

unset(_amrex_cuda_archs)
unset(_amrex_cuda_archs_legacy)
unset(_amrex_cuda_archs_given)
unset(_amrex_cuda_archs_cmdline)
unset(_amrex_cuda_archs_cmdline_now)
unset(_amrex_cuda_archs_cached)
unset(_amrex_cuda_archs_had_cache)
unset(_amrex_cuda_archs_shadowed)
unset(_amrex_cuda_arch_cmdline_now)
