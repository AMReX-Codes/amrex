if (NOT ENABLE_FORTRAN_INTERFACES)
    return()
endif ()

if (DIM EQUAL 1)
    return()
endif ()

#
# Use directory name as exe name
#
get_filename_component(exe_name ${CMAKE_CURRENT_LIST_DIR} NAME)
set(exe_name "${exe_name}.exe")

#
# This test works for both 2D and 3D builds
#
add_executable( ${exe_name} EXCLUDE_FROM_ALL "")

target_sources( ${exe_name}
    PRIVATE
    main.F90
    mytest.F90
    init_prob.F90
    )

set_target_properties( ${exe_name}
    PROPERTIES
    INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
    Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
    RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
    )

target_link_libraries( ${exe_name} PRIVATE amrex )

# Find input files
file( GLOB_RECURSE inputs LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input* )

# Copy input files to corresponding build dir
file( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# Add to the "tutorial" target
add_tutorial(${exe_name})
