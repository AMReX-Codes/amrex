add_sources( AMReX_AmrData.H        AMReX_XYPlotDataList.H )
add_sources( AMReX_DataServices.H   AMReX_AmrvisConstants.H )
add_sources( AMReX_AmrData.cpp      AMReX_XYPlotDataList.cpp )
add_sources( AMReX_DataServices.cpp AMReX_FABUTIL_${DIM}D.F )

if (ENABLE_PROPARSER)
   add_sources( AMReX_XYPlotDataList.cpp AMReX_XYPlotDataList.H )
endif ()
