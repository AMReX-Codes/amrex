set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

set ( CXXSRC
   AMReX_BndryData.cpp  AMReX_BndryRegister.cpp  AMReX_FabSet.cpp
   AMReX_InterpBndryData.cpp  AMReX_MacBndry.cpp  AMReX_Mask.cpp
   AMReX_MultiMask.cpp )

set ( F77SRC
   AMReX_INTERPBNDRYDATA_${BL_SPACEDIM}D.F  AMReX_LO_UTIL.F )

set ( ALLHEADERS
   AMReX_BndryData.H  AMReX_BndryRegister.H  AMReX_BoundCond.H
   AMReX_FabSet.H  AMReX_INTERPBNDRYDATA_F.H  AMReX_InterpBndryData.H
   AMReX_LO_BCTYPES.H  AMReX_MacBndry.H  AMReX_Mask.H  AMReX_MultiMask.H )


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
