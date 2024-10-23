add_amrex_define(AMREX_USE_FFT NO_LEGACY)

foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_FFT.H
       AMReX_FFT.cpp
       AMReX_FFT_R2C.H
       AMReX_FFT_R2X.H
       AMReX_FFT_Helper.H
       AMReX_FFT_Poisson.H
       )

endforeach()
