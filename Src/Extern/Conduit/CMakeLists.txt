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
set_target_properties(conduit:conduit PROPERTIES IMPORTED_GLOBAL True )

if (ENABLE_MPI)
   target_link_libraries(amrex PUBLIC conduit:conduit_mpi)
   set_target_properties(conduit:conduit_mpi PROPERTIES IMPORTED_GLOBAL True )   
endif ()

if (ENABLE_ASCENT)
   find_package(Ascent)
   
   if (ENABLE_MPI)
      target_link_libraries( amrex PUBLIC ascent::ascent_mpi )
      set_target_properties( ascent::ascent_mpi PROPERTIES IMPORTED_GLOBAL True )
   else ()
      target_link_libraries( amrex PUBLIC ascent::ascent )
      set_target_properties( ascent::ascent PROPERTIES IMPORTED_GLOBAL True )
   endif ()   
endif ()



