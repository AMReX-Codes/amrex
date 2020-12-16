#
# HDF5 -- here it would be best to create an imported target
#
if (AMReX_HDF5)
    set(HDF5_PREFER_PARALLEL TRUE)
    find_package(HDF5 1.10.4 REQUIRED)
    if (NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR "\nHDF5 library does not support parallel I/O")
     endif ()

    if (TARGET hdf5::hdf5)  # CMake >= 3.19
       target_link_libraries(amrex PUBLIC hdf5::hdf5)
    else ()  # CMake < 3.19 -- Remove when minimum cmake version is bumped up
       target_include_directories(amrex PUBLIC ${HDF5_INCLUDE_DIRS})
       target_compile_definitions(amrex PUBLIC ${HDF5_DEFINITIONS})
       target_link_libraries(amrex PUBLIC ${HDF5_LIBRARIES})
    endif ()

endif ()


#
# Sensei
#
if (AMReX_SENSEI)
    find_package(SENSEI REQUIRED)
    target_link_libraries( amrex PUBLIC sensei )
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
    find_package(HYPRE 2.18.2 REQUIRED)
    target_link_libraries( amrex PUBLIC HYPRE )
endif ()


#
# PETSc
#
if (AMReX_PETSC)
    find_package(PETSc 2.13 REQUIRED)
    target_link_libraries( amrex PUBLIC PETSC )
endif ()
