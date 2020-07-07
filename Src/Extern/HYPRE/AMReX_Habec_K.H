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

AMREX_FORCE_INLINE
void amrex_hpdiag (Box const& box,
                   BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                   FArrayBox& diag)
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real> const& diag_arr = diag.array();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        diag_arr(i,j,k) = 1.0/mat_arr(i,j,k)[0];
        for (int ii=0; ii<2*AMREX_SPACEDIM + 1; ii++){
            mat_arr(i,j,k)[ii] = mat_arr(i,j,k)[ii]*diag_arr(i,j,k);
        }
    });
}

#ifdef AMREX_USE_EB

template<typename T> AMREX_FORCE_INLINE
void amrex_hpeb_fill_cellid (Box const& box,
                             T& nrows,
                             BaseFab<T>& cell_id,
                             const EBCellFlagFab& flag)
{
    Array4<T> const& cell_id_arr = cell_id.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            cell_id_arr(i,j,k) = *nrowsg;
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
    nrows = nrows_gpu.dataValue();
}

AMREX_FORCE_INLINE
void amrex_hpeb_copy_from_vec (Box const& box,
                               FArrayBox& a,
                               Real* v,
                               const EBCellFlagFab& flag)
{
    Array4<Real> const& a_arr = a.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    int nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            a_arr(i,j,k) = v[*nrowsg];
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
}

AMREX_FORCE_INLINE
void amrex_hpeb_copy_to_vec (Box const& box,
                             FArrayBox& a,
                             Real* v,
                             const EBCellFlagFab& flag)
{
    Array4<Real> const& a_arr = a.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    int nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            v[*nrowsg] = a_arr(i,j,k);
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
}

#endif

}
#endif
