target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# TODO: Particles PR merges another file
target_sources(amrex
   PRIVATE 
   AMReX_Conduit_Blueprint.H
   AMReX_Conduit_Blueprint.cpp
   )

#
# Dependencies
#

# For now do not look for specific components (if Conduit has any)
# It would be good to check whether component conduit_mpi is available
# or not
find_package(Conduit REQUIRED)

target_link_libraries(amrex PUBLIC conduit:conduit)

if (ENABLE_MPI)
   target_link_libraries(amrex PUBLIC conduit:conduit_mpi)   
endif ()

if (ENABLE_ASCENT)
   find_package(Ascent)
   
   if (ENABLE_MPI)
      target_link_libraries( amrex PUBLIC ascent::ascent_mpi )
   else ()
      target_link_libraries( amrex PUBLIC ascent::ascent )
   endif ()
   
endif ()

