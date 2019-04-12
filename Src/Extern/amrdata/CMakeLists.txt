target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE
   AMReX_AmrData.H
   AMReX_AmrData.cpp
   AMReX_XYPlotDataList.H
   AMReX_XYPlotDataList.cpp
   AMReX_DataServices.H
   AMReX_DataServices.cpp 
   AMReX_AmrvisConstants.H
   AMReX_FABUTIL_${DIM}D.F
   )

if (ENABLE_PROPARSER)
   target_sources( amrex PRIVATE AMReX_XYPlotDataList.cpp AMReX_XYPlotDataList.H )
endif ()
