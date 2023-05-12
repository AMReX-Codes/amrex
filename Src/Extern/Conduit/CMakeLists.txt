foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}D PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    # TODO: Particles PR merges another file
    target_sources(amrex_${D}D
       PRIVATE
       AMReX_Conduit_Blueprint.H
       AMReX_Conduit_Blueprint.cpp
    )
endforeach()
