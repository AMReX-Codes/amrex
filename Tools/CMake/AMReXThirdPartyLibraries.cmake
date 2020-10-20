#
# HDF5 -- here it would be best to create an imported target
#
if (AMReX_HDF5)
    set(HDF5_PREFER_PARALLEL TRUE)
    find_package(HDF5 1.10.4 REQUIRED COMPONENTS CXX)
    if (NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR "\nHDF5 library does not support parallel I/O")
    endif ()
    target_include_directories(amrex PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
    target_compile_definitions(amrex PUBLIC ${HDF5_CXX_DEFINES})
    target_link_libraries(amrex PUBLIC ${HDF5_CXX_LIBRARIES})
endif ()


#
# Sensei
#
if (AMReX_SENSEI)
    find_package(SENSEI REQUIRED)
    target_link_libraries( amrex PUBLIC sensei )
endif ()


#
# SUNDIALS
#
if (AMReX_SUNDIALS)
    # We link to libraries and always include nvecserial (in case app code needs it)
    set(_sundials_components farkode_mod;fcvode_mod)

    find_package(SUNDIALS 4 REQUIRED COMPONENTS ${_sundials_components})

    foreach (_comp ${_sundials_components})
        target_link_libraries(amrex PUBLIC SUNDIALS::${_comp})
    endforeach ()

    unset(_sundials_components)
endif ()


#
#  Ascent
#
if (AMReX_ASCENT) # Ascent will find conduit, so check for Ascent first
    find_package(Ascent REQUIRED)
    if (AMReX_MPI)
        target_link_libraries( amrex PUBLIC ascent::ascent_mpi )
    else ()
        target_link_libraries( amrex PUBLIC ascent::ascent )
    endif ()
endif ()


#
# Conduit
#
if (AMReX_CONDUIT)
    find_package(Conduit REQUIRED)
    if (AMReX_MPI)
        target_link_libraries( amrex PUBLIC conduit::conduit_mpi )
    else ()
        target_link_libraries( amrex PUBLIC conduit::conduit )
    endif ()
endif ()


#
# HYPRE
#
if (AMReX_HYPRE)
    find_package(HYPRE 2.15 REQUIRED)
    target_link_libraries( amrex PUBLIC HYPRE )
endif ()


#
# PETSc
#
if (AMReX_PETSC)
    find_package(PETSc 2.13 REQUIRED)
    target_link_libraries( amrex PUBLIC PETSC )
endif ()
