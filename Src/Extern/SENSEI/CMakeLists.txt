add_sources ( AMReX_AmrMeshDataAdaptor.H AMReX_AmrDataAdaptor.H
    AMReX_InSituBridge.H AMReX_AmrMeshInSituBridge.H AMReX_AmrInSituBridge.H
    AMReX_InSituUtils.H )

add_sources ( AMReX_AmrMeshDataAdaptor.cpp AMReX_AmrDataAdaptor.cpp
    AMReX_InSituBridge.cpp AMReX_AmrMeshInSituBridge.cpp AMReX_AmrInSituBridge.cpp
    AMReX_InSituUtils.cpp )

target_link_libraries( amrex PUBLIC sensei )
