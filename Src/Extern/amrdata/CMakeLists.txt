foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_AmrData.H
       AMReX_AmrData.cpp
       AMReX_XYPlotDataList.H
       AMReX_XYPlotDataList.cpp
       AMReX_DataServices.H
       AMReX_DataServices.cpp
       AMReX_WritePlotFile.H
       AMReX_WritePlotFile.cpp
       AMReX_AmrvisConstants.H
       )

    if (AMReX_PROFPARSER)
       target_sources(amrex_${D}d PRIVATE AMReX_XYPlotDataList.cpp AMReX_XYPlotDataList.H )
    endif ()
endforeach()
