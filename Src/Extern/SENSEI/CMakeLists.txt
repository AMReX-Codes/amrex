#
# Compile definitions required by Extern/SENSEI
# 
add_amrex_define(BL_USE_SENSEI_INSITU)

#
# Source files
# 
add_sources( AMReX_AmrMeshDataAdaptor.H AMReX_AmrDataAdaptor.H
    AMReX_InSituBridge.H AMReX_AmrMeshInSituBridge.H AMReX_AmrInSituBridge.H
    AMReX_InSituUtils.H )

add_sources( AMReX_AmrMeshDataAdaptor.cpp AMReX_AmrDataAdaptor.cpp
    AMReX_InSituBridge.cpp AMReX_AmrMeshInSituBridge.cpp AMReX_AmrInSituBridge.cpp
    AMReX_InSituUtils.cpp )

# This is not gonna work since "sensei" is not a target of this project
# If it is a target in some other project, said project should be properly found/imported
# in CMake
target_link_libraries( amrex PUBLIC sensei )
