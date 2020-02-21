target_compile_definitions(amrex PUBLIC $<BUILD_INTERFACE:AMREX_USE_PROFPARSER BL_USE_PROFPARSER>)
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
    PRIVATE
    AMReX_BLProfStats.H
    AMReX_BLProfStats.cpp
    AMReX_BLProfUtilities.H
    AMReX_BLProfUtilities.cpp
    AMReX_CommProfStats.cpp
    AMReX_RegionsProfStats.cpp
    AMReX_AVGDOWN_${DIM}D.F
    )
