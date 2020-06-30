#ifndef AMREX_Habec_K_H_
#define AMREX_Habec_K_H_

#if (AMREX_SPACEDIM == 2)
#include <AMReX_Habec_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_Habec_3D_K.H>
#endif

namespace amrex {

AMREX_FORCE_INLINE
void amrex_hpacoef (Box const& box,
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& a,
                    Real& sa)
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& a_arr = a.const_array();
    if (sa == 0.0){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = 0.0;
        });
    }else{
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = sa * a_arr(i,j,k);
        });
    }
}

}
#endif
