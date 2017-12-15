set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND CXXSRC     AMReX_FabSet.cpp AMReX_BndryRegister.cpp AMReX_Mask.cpp AMReX_MultiMask.cpp )
list ( APPEND CXXSRC     AMReX_BndryData.cpp AMReX_InterpBndryData.cpp AMReX_MacBndry.cpp )

list ( APPEND ALLHEADERS AMReX_FabSet.H AMReX_BndryRegister.H AMReX_Mask.H AMReX_MultiMask.H )
list ( APPEND ALLHEADERS AMReX_BndryData.H AMReX_BoundCond.H AMReX_InterpBndryData.H )
list ( APPEND ALLHEADERS AMReX_LO_BCTYPES.H AMReX_MacBndry.H )

list ( APPEND ALLHEADERS AMReX_INTERPBNDRYDATA_F.H )
list ( APPEND F77SRC     AMReX_INTERPBNDRYDATA_${DIM}D.F AMReX_LO_UTIL.F )

list ( APPEND ALLHEADERS AMReX_YAFluxRegister.H   AMReX_YAFluxRegister_F.H )
list ( APPEND CXXSRC     AMReX_YAFluxRegister.cpp )
list ( APPEND F90SRC     AMReX_YAFluxRegister_${DIM}d.F90 AMReX_YAFluxRegister_nd.F90 )

list ( APPEND F90SRC     AMReX_lo_bctypes_mod.F90 )
# 
# Accumulate sources
# 
set ( ALLSRC ${CXXSRC} ${F77SRC} ${F90SRC} )

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

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
