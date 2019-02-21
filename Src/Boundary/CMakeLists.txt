#
# Add all sources
# 
add_sources ( AMReX_FabSet.cpp AMReX_BndryRegister.cpp AMReX_Mask.cpp AMReX_MultiMask.cpp )
add_sources ( AMReX_BndryData.cpp AMReX_InterpBndryData.cpp AMReX_MacBndry.cpp )

add_sources ( AMReX_FabSet.H AMReX_BndryRegister.H AMReX_Mask.H AMReX_MultiMask.H )
add_sources ( AMReX_BndryData.H AMReX_BoundCond.H AMReX_InterpBndryData.H )
add_sources ( AMReX_LO_BCTYPES.H AMReX_MacBndry.H )

add_sources ( AMReX_InterpBndryData_K.H AMReX_InterpBndryData_${DIM}D_K.H )

add_sources ( AMReX_LO_UTIL.F90 )
add_sources ( AMReX_LOUtil_K.H )

add_sources ( AMReX_YAFluxRegister.H )
add_sources ( AMReX_YAFluxRegister_K.H AMReX_YAFluxRegister_${DIM}D_K.H )
add_sources ( AMReX_YAFluxRegister.cpp )

add_sources ( AMReX_lo_bctypes_mod.F90 )

