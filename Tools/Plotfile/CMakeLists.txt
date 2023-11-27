project(AMReX-PlotfileTools
   DESCRIPTION "AMReX tools to postprocess plot files"
   )

# Plotfile tools requires only AMReX sources from Src/Base which get always
# compiled, so no need to check which part of AMReX is enabled by user

# List of plotfile targets
set(_exe_names
   fboxinfo
   fcompare
   fextract
   fextrema
   fgradient
   fnan
   fsnapshot
   ftime
   fvarnames
   fvolumesum
   )

# Build targets one by one
foreach( _exe IN LISTS _exe_names)
   add_executable(${_exe} ${_exe}.cpp)
   add_executable(AMReX::${_exe} ALIAS ${_exe})
   target_link_libraries(${_exe} PRIVATE amrex)
   if (AMReX_CUDA)
      set_source_files_properties(${_exe}.cpp PROPERTIES LANGUAGE CUDA)
   endif()

    # Add prefix to each tool's name to make them unique when installed.
    # This avoids potential collisions of names on user systems, e.g., in
    # software packages (Spack/Conda/Debian/...).
    set_target_properties(${_exe} PROPERTIES OUTPUT_NAME "amrex_${_exe}")
endforeach()


# target snapshot needs a special treatment
target_include_directories(fsnapshot PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(fsnapshot PRIVATE AMReX_PPMUtil.H AMReX_PPMUtil.cpp)
if (AMReX_CUDA)
   set_source_files_properties(AMReX_PPMUtil.cpp PROPERTIES LANGUAGE CUDA)
   target_compile_features(fsnapshot PUBLIC cxx_std_17)
endif()
