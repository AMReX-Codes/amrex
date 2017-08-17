set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

# 
# For now we accumulate sources for folder amrdata only
#

set (CXXSRC amrdata/AMReX_AmrData.cpp  amrdata/AMReX_DataServices.cpp)

set (F77SRC amrdata/AMReX_FABUTIL_${BL_SPACEDIM}D.F)

set (ALLHEADERS
   amrdata/AMReX_AmrData.H  amrdata/AMReX_AmrvisConstants.H  amrdata/AMReX_DataServices.H )


# Accumulate sources
set ( ALLSRC ${CXXSRC} ${F77SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

