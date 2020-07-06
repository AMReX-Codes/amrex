#ifndef AMREX_HABEC_2D_H_
#define AMREX_HABEC_2D_H_

#include <AMReX_LO_BCTYPES.H>
#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBABecLap_K.H>
#endif

namespace amrex {

AMREX_FORCE_INLINE
void amrex_hpbcoef (Box const& box,
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& b,
                    Real& sb,
                    const Real* dx,
                    int& idim)
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& b_arr = b.const_array();
    Real fac = sb / (dx[idim]*dx[idim]);

    if (idim == 0) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + fac * (b_arr(i,j,k) + b_arr(i+1,j,k));
            mat_arr(i,j,k)[1] = - fac * b_arr(i,j,k);
            mat_arr(i,j,k)[2] = - fac * b_arr(i+1,j,k);
        });
    }else{
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + fac * (b_arr(i,j,k) + b_arr(i,j+1,k));
            mat_arr(i,j,k)[3] = - fac * b_arr(i,j,k);
            mat_arr(i,j,k)[4] = - fac * b_arr(i,j+1,k);
        });
    }
}

AMREX_FORCE_INLINE
void amrex_hpmat (Box const& box,
                  BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                  const FArrayBox& b,
                  const Mask& msk,
                  Real& sb,
                  const Real* dx,
                  int& cdir,
                  const int& bct,
                  const Real& bcl,
                  const int& bho)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& b_arr = b.const_array();
    Array4<int const> const& msk_arr = msk.const_array();

    Real h = 0.0;
    if (cdir == 0 || cdir == 2) {
        h = dx[0];
    }else{
        h = dx[1];
    }
    Real fac = sb/(h*h);

    Real h2 = 0.0;
    Real h3 = 0.0;
    Real bf1 = 0.0;
    Real bf2 = 0.0;
    if (bct == AMREX_LO_DIRICHLET) {
        h2 = 0.5*h;
        if (bho >= 1) {
            h3 = 3.0*h2;
            bf1 = fac*((h3 - bcl)/(bcl + h2) - 1.0);
            bf2 = fac*(bcl - h2)/(bcl + h3);
        }else{
            bf1 = fac*( h/(bcl + h2) - 1.0);
            bf2 = 0.0;
        }
    }else if (bct == AMREX_LO_NEUMANN) {
        bf1 = -fac;
        bf2 = 0.0;
    }else{
        Abort("hpmat: unsupported boundary type");
    }

    if (cdir == 0){
        Box tbx = box;
        tbx.setSmall(0,lo.x).setBig(0,lo.x); // 0: x-direction
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (tbx, i, j, k,
        {
            if (msk_arr(lo.x-1,j,k) > 0){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k);
                mat_arr(i,j,k)[1] = 0.0;
                mat_arr(i,j,k)[2] = mat_arr(i,j,k)[2] + bf2*b_arr(i,j,k);
            }
        });
    }else if (cdir == 2){
        Box tbx = box;
        tbx.setSmall(0,hi.x).setBig(0,hi.x); // 0: x-direction
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (tbx, i, j, k,
        {
            if (msk_arr(hi.x+1,j,k) > 0){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i+1,j,k);
                mat_arr(i,j,k)[2] = 0.0;
                mat_arr(i,j,k)[1] = mat_arr(i,j,k)[1] + bf2*b_arr(i+1,j,k);
            }
        });
   }else if (cdir == 1){
        Box tby = box;
        tby.setSmall(1,lo.y).setBig(1,lo.y); // 1: y-direction
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (tby, i, j, k,
        {
            if (msk_arr(i,lo.y-1,k) > 0){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k);
                mat_arr(i,j,k)[3] = 0.0;
                mat_arr(i,j,k)[4] = mat_arr(i,j,k)[4] + bf2*b_arr(i,j,k);
            }
        });
   }else if (cdir == 3){
        Box tby = box;
        tby.setSmall(1,hi.y).setBig(1,hi.y); // 1: y-direction
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (tby, i, j, k,
        {
            if (msk_arr(i,hi.y+1,k) > 0){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j+1,k);
                mat_arr(i,j,k)[4] = 0.0;
                mat_arr(i,j,k)[3] = mat_arr(i,j,k)[3] + bf2*b_arr(i,j+1,k);
            }
        });
   }else{
        Abort("hpmat: impossible face orientation");
   }
}

template<typename T> AMREX_FORCE_INLINE
void amrex_hpijmatrix (Box const& box,
                       const T& nrows, 
                       T* ncols,
                       T* rows, 
                       T* colsg,
                       Real* matg, 
                       BaseFab<T>& cell_id,
                       T& offset, FArrayBox& diaginv,
                       FArrayBox& acoefs, FArrayBox& bcoefsx,
                       FArrayBox& bcoefsy,
                       Real& sa, Real& sb, const Real* dx, 
                       GpuArray<int,AMREX_SPACEDIM*2>& bct,
                       GpuArray<Real,AMREX_SPACEDIM*2> bcl,
                       const int& bho)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    GpuArray<Real,AMREX_SPACEDIM> fac;
    GpuArray<Real,AMREX_SPACEDIM*2> bf1;
    GpuArray<Real,AMREX_SPACEDIM*2> bf2;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        fac[i] = sb/(dx[i]*dx[i]);
    }
   
    int idim = 0;
    for (int cdir=0; cdir<=AMREX_SPACEDIM*2-1; cdir++){ 
       if (cdir == 0 || cdir == 2){
          idim = 0;
       }else{
          idim = 1;
       }

       Real h = dx[idim];

       if (bct[cdir] == AMREX_LO_DIRICHLET) {
          Real h2 = 0.5*h;
          if (bho>=1) {
             Real h3 = 3.0*h2;
             bf1[cdir] = fac[idim] * ((h3 - bcl[cdir]) / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = fac[idim] * (bcl[cdir] - h2) / (bcl[cdir] + h3);
          }else{
             bf1[cdir] = fac[idim] * ( h / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = 0.0;
          }
       }else if (bct[cdir] == AMREX_LO_NEUMANN) {
          bf1[cdir] = -fac[idim];
          bf2[cdir] = 0.0;
       }
    }

    Gpu::DeviceVector<int> iter(2,0);

    Array4<T> const& cell_id_arr = cell_id.array();
    Gpu::DeviceVector<T> cols_tmp(AMREX_SPACEDIM*2+1,0.0);
    Gpu::DeviceVector<Real> mat_tmp(AMREX_SPACEDIM*2+1,0.0);
    Array4<Real const> const& a_arr = acoefs.const_array();
    Array4<Real const> const& bx_arr = bcoefsx.const_array();
    Array4<Real const> const& by_arr = bcoefsy.const_array();
    Array4<Real> const& diag_arr = diaginv.array();

    auto cols_tmpPtr = cols_tmp.dataPtr();
    auto mat_tmpPtr = mat_tmp.dataPtr();
    auto iterPtr = iter.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        int irow = iterPtr[0];
        rows[irow]  = cell_id_arr(i,j,k);
        ncols[irow] = 0;

        cols_tmpPtr[0] = cell_id_arr(i,j,k);
        mat_tmpPtr[0]  = sa*a_arr(i,j,k) + fac[0]*(bx_arr(i,j,k)+bx_arr(i+1,j,k)) 
                                         + fac[1]*(by_arr(i,j,k)+by_arr(i,j+1,k));

        cols_tmpPtr[1] = cell_id_arr(i-1,j,k);
        mat_tmpPtr[1]  = -fac[0]*bx_arr(i,j,k);

        cols_tmpPtr[2] = cell_id_arr(i+1,j,k);
        mat_tmpPtr[2]  = -fac[0]*bx_arr(i+1,j,k);

        cols_tmpPtr[3] = cell_id_arr(i,j-1,k);
        mat_tmpPtr[3]  = -fac[1]*by_arr(i,j,k);

        cols_tmpPtr[4] = cell_id_arr(i,j+1,k);
        mat_tmpPtr[4]  = -fac[1]*by_arr(i,j+1,k);

        if (i == lo.x && cell_id_arr(i-1,j,k)<0) {
           mat_tmpPtr[0] += bf1[0]*bx_arr(i,j,k);
           mat_tmpPtr[2] += bf2[0]*bx_arr(i,j,k);
        }

        if (i == hi.x && cell_id_arr(i+1,j,k)<0) {
           mat_tmpPtr[0] += bf1[2]*bx_arr(i+1,j,k);
           mat_tmpPtr[1] += bf2[2]*bx_arr(i+1,j,k);
        }

        if (j == lo.y && cell_id_arr(i,j-1,k)<0) {
           mat_tmpPtr[0] += bf1[1]*by_arr(i,j,k);
           mat_tmpPtr[4] += bf2[1]*by_arr(i,j,k);
        }

        if (j == hi.y && cell_id_arr(i,j+1,k)<0) {
           mat_tmpPtr[0] += bf1[3]*by_arr(i,j+1,k);
           mat_tmpPtr[3] += bf2[3]*by_arr(i,j+1,k);
        }

        diag_arr(i,j,k) = 1.0/mat_tmpPtr[0];

        for (int ic=0; ic<=4; ic++){
            if (cols_tmpPtr[ic] >= 0){
                int imat = iterPtr[1];
                ncols[irow]  = ncols[irow] +1;
                colsg[imat]  = cols_tmpPtr[ic];
                matg[imat]   = mat_tmpPtr[ic]*diag_arr(i,j,k);
                Gpu::Atomic::Add(&iterPtr[1], 1);
            }
        }
        Gpu::Atomic::Add(&iterPtr[0], 1);
    });
}

#ifdef AMREX_USE_EB

template<typename T> AMREX_FORCE_INLINE
void amrex_hpeb_ijmatrix (Box const& box,
                          const T& nrows,
                          T* ncols,
                          T* rows,
                          T* colsg,
                          Real* matg,
                          BaseFab<T>& cell_id,
                          T& offset, FArrayBox& diaginv,
                          FArrayBox& acoefs, FArrayBox& bcoefsx,
                          FArrayBox& bcoefsy,
                          const EBCellFlagFab& flag,
                          const FArrayBox& vfrc,
                          const FArrayBox& apx, const FArrayBox& apy,
                          const FArrayBox& fcx,
                          const FArrayBox& fcy,
                          const FArrayBox& ba, const FArrayBox& bcen,
                          const FArrayBox& beb,
                          const int& is_eb_dirichlet,
                          Real& sa, Real& sb, const Real* dx,
                          GpuArray<int,AMREX_SPACEDIM*2>& bct,
                          GpuArray<Real,AMREX_SPACEDIM*2> bcl,
                          const int& bho)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    GpuArray<Real,AMREX_SPACEDIM> fac;
    GpuArray<Real,AMREX_SPACEDIM*2> bf1;
    GpuArray<Real,AMREX_SPACEDIM*2> bf2;
    GpuArray<Real,AMREX_SPACEDIM*2> bflo;
 
    int is_dirichlet = (is_eb_dirichlet != 0) ? 1 : 0;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        fac[i] = sb/(dx[i]*dx[i]);
    }
   
    int idim = 0;
    for (int cdir=0; cdir<2*AMREX_SPACEDIM; cdir++){ 
       if (cdir == 0 || cdir == 2){
          idim = 0;
       }else{
          idim = 1;
       }

       Real h = dx[idim];

       if (bct[cdir] == AMREX_LO_DIRICHLET) {
          Real h2 = 0.5*h;
          bflo[cdir] = fac[idim] * ( h / (bcl[cdir] + h2) - 1.0); 
          if (bho>=1) {
             Real h3 = 3.0*h2;
             bf1[cdir] = fac[idim] * ((h3 - bcl[cdir]) / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = fac[idim] * (bcl[cdir] - h2) / (bcl[cdir] + h3);
          }else{
             bf1[cdir] = bflo[cdir];
             bf2[cdir] = 0.0;
          }
       }else if (bct[cdir] == AMREX_LO_NEUMANN) {
          bflo[cdir] = -fac[idim];
          bf1[cdir] = -fac[idim];
          bf2[cdir] = 0.0;
       }
    }

    Gpu::DeviceVector<int> iter(2,0);
    Array4<T> const& cell_id_arr = cell_id.array();
    Array4<Real const> const& a_arr = acoefs.const_array();
    Array4<Real const> const& bx_arr = bcoefsx.const_array();
    Array4<Real const> const& by_arr = bcoefsy.const_array();
    Array4<Real> const& diag_arr = diaginv.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();
    Array4<Real const> const& apx_arr = apx.const_array();
    Array4<Real const> const& fcx_arr = fcx.const_array();
    Array4<Real const> const& apy_arr = apy.const_array();
    Array4<Real const> const& fcy_arr = fcy.const_array();
    Array4<Real const> const& bcen_arr = bcen.const_array();
    Array4<Real const> const& vfrc_arr = vfrc.const_array();
    Array4<Real const> const& ba_arr = ba.const_array();
    Array4<Real const> const& beb_arr = beb.const_array();

    auto iterPtr = iter.dataPtr();

    AMREX_HOST_DEVICE_FOR_3D (box, i, j, k,
    {
        hpeb_ij(i,j,k,nrows,ncols,rows,colsg,matg,fac,bf1,bf2,bflo,
                  lo,hi,sa,bho,is_dirichlet,
                  iterPtr,cell_id_arr,a_arr,bx_arr,
                  by_arr,diag_arr,flag_arr,apx_arr,fcx_arr,
                  apy_arr,fcy_arr,bcen_arr,vfrc_arr,
                  ba_arr,beb_arr);
    });

}

template<typename T> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hpeb_ij (int i, int j, int k,
              const T& nrows,
              T* ncols,
              T* rows,
              T* colsg,
              Real* matg,
              const GpuArray<Real,AMREX_SPACEDIM>& fac,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bf1,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bf2,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bflo,
              const Dim3 lo,
              const Dim3 hi,
              const Real& sa, const int& bho, const int& is_dirichlet,
              int* iterPtr,
              Array4<T> const& cell_id_arr,
              Array4<Real const> const& a_arr,
              Array4<Real const> const& bx_arr,
              Array4<Real const> const& by_arr,
              Array4<Real> const& diag_arr,
              Array4<const EBCellFlag> const& flag_arr,
              Array4<Real const> const& apx_arr,
              Array4<Real const> const& fcx_arr,
              Array4<Real const> const& apy_arr,
              Array4<Real const> const& fcy_arr,
              Array4<Real const> const& bcen_arr,
              Array4<Real const> const& vfrc_arr,
              Array4<Real const> const& ba_arr,
              Array4<Real const> const& beb_arr) noexcept
{
        Array2D<Real,-1,1,-1,1> mat_tmp;
        Array1D<Real,0,3> phig1;
        Array1D<Real,0,3> feb;

        if (flag_arr(i,j,k).isCovered()){
            diag_arr(i,j,k) = 0.0;
        }else{
            int irow = iterPtr[0];
            rows[irow]  = cell_id_arr(i,j,k);
            ncols[irow] = 0;
            for (int jj=-1; jj<=1; jj++){
                for (int ii=-1; ii<=1; ii++){
                     mat_tmp(ii,jj) = 0.0;
                }
            }

            if (flag_arr(i,j,k).isRegular()){
                mat_tmp(0,0) = sa*a_arr(i,j,k) + fac[0]*(bx_arr(i,j,k)+bx_arr(i+1,j,k))
                                                 + fac[1]*(by_arr(i,j,k)+by_arr(i,j+1,k));
                mat_tmp(-1, 0) = -fac[0]*bx_arr(i,j,k);
                mat_tmp( 1, 0) = -fac[0]*bx_arr(i+1,j,k);
                mat_tmp( 0,-1) = -fac[1]*by_arr(i,j,k);
                mat_tmp( 0, 1) = -fac[1]*by_arr(i,j+1,k);

                if (i == lo.x && cell_id_arr(i-1,j,k)<0) {
                   mat_tmp(0,0)  = mat_tmp(0,0) + bf1[0]*bx_arr(i,j,k);
                   mat_tmp(-1,0) = 0.0;
                   mat_tmp(1,0)  = mat_tmp(1,0) + bf2[0]*bx_arr(i,j,k);
                }

                if (i == hi.x && cell_id_arr(i+1,j,k)<0) {
                   mat_tmp(0,0)  = mat_tmp(0,0) + bf1[2]*bx_arr(i+1,j,k);
                   mat_tmp(1,0)  = 0.0;
                   mat_tmp(-1,0) = mat_tmp(-1,0) + bf2[2]*bx_arr(i+1,j,k);
                }

                if (j == lo.y && cell_id_arr(i,j-1,k)<0) {
                   mat_tmp(0,0)  = mat_tmp(0,0) + bf1[1]*by_arr(i,j,k);
                   mat_tmp(0,-1) = 0.0;
                   mat_tmp(0,1)  = mat_tmp(0,1) + bf2[1]*by_arr(i,j,k);
                }

                if (j == hi.y && cell_id_arr(i,j+1,k)<0) {
                   mat_tmp(0,0) = mat_tmp(0,0) + bf1[3]*by_arr(i,j+1,k);
                   mat_tmp(0,1) = 0.0;
                   mat_tmp(0,-1) = mat_tmp(0,-1) + bf2[3]*by_arr(i,j+1,k);
                }

            }else{
                   int cdir = 0;
                   int joff = 0;
                   int koff = 0;
                   int jj = 0;
                   int kk = 0;
                   Real f = fac[0];
                   Real area = apx_arr(i,j,k);
                   Real bc = bx_arr(i,j,k);
                   Real fracy = 0.0;
                   if (area > 0.0) {
                      if (area != 1.0) {
                         joff = static_cast<int>(Math::copysign(1.0, fcx_arr(i,j,k)));
                         jj = j+joff;
                         if (cell_id_arr(i-1,jj,k)<0 && cell_id_arr(i,jj,k)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcx_arr(i,j,k));
                         }
                      }else{
                         joff = 0;
                         jj = j;
                         fracy = 0.0;
                      }

                      Real tmp = (1.0-fracy)*area*bc;
                      if (cell_id_arr(i-1,j,k)>=0) {
                         mat_tmp( 0,0) = mat_tmp( 0,0) + tmp*f;
                         mat_tmp(-1,0) = mat_tmp(-1,0) - tmp*f;
                      }else if (cell_id_arr(i+1,j,k)<0 || apx_arr(i+1,j,k)==0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(1,0) = mat_tmp(1,0) + tmp*   bf2[cdir];
                      }

                      if (fracy>0.0) {
                         if (cell_id_arr(i-1,jj,k)>=0 && cell_id_arr(i,jj,k)>=0) {
                            mat_tmp(-1,joff) = mat_tmp(-1,joff) - fracy*area*bx_arr(i,jj,k)*f;
                            mat_tmp( 0,joff) = mat_tmp( 0,joff) + fracy*area*bx_arr(i,jj,k)*f;
                         }else if (cell_id_arr(i+1,jj,k)<0 || apx_arr(i+1,jj,k)==0) {
                            mat_tmp( 0,joff) = mat_tmp( 0,joff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff) = mat_tmp( 0,joff) + tmp*(f+bf1[cdir]);
                            mat_tmp( 1,joff) = mat_tmp( 1,joff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 2;
                   f = fac[0];
                   area = apx_arr(i+1,j,k);
                   bc = bx_arr(i+1,j,k);
                   if (area>0.0) {
                      if (area!=1.0) {
                         joff = static_cast<int>(Math::copysign(1.0, fcx_arr(i+1,j,k)));
                         jj = j+joff;
                         if (cell_id_arr(i,jj,k)<0 && cell_id_arr(i+1,jj,k)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcx_arr(i+1,j,k));
                         }
                      }else{
                         joff = 0;
                         jj = j;
                         fracy = 0.0;
                      }

                      Real tmp = (1.0-fracy)*area*bc;
                      if (cell_id_arr(i+1,j,k)>=0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*f;
                         mat_tmp(1,0) = mat_tmp(1,0) - tmp*f;
                      }else if (cell_id_arr(i-1,j,k)<0 || apx_arr(i,j,k)==0.0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp( 0,0) = mat_tmp( 0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(-1,0) = mat_tmp(-1,0) + tmp*   bf2[cdir];
                      }

                      if (fracy>0.0) {
                         if (cell_id_arr(i,jj,k)>=0 && cell_id_arr(i+1,jj,k)>=0) {
                            mat_tmp(0,joff) = mat_tmp(0,joff) + fracy*area*bx_arr(i+1,jj,k)*f;
                            mat_tmp(1,joff) = mat_tmp(1,joff) - fracy*area*bx_arr(i+1,jj,k)*f;
                         }else if (cell_id_arr(i-1,jj,k)<0 || apx_arr(i,jj,k)==0.0 || bho==0) {
                            mat_tmp(0,joff) = mat_tmp(0,joff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff) = mat_tmp( 0,joff) + tmp*(f+bf1[cdir]);
                            mat_tmp(-1,joff) = mat_tmp(-1,joff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 1;
                   Real fracx = 0.0;
                   int ii = 0;
                   int ioff = 0;

                   f = fac[1];
                   area = apy_arr(i,j,k);
                   bc = by_arr(i,j,k);
                   if (area>0) {
                      if (area!=1) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j,k)));
                         ii = i+ioff;
                         if (cell_id_arr(ii,j-1,k)<0 && cell_id_arr(ii,j,k)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcy_arr(i,j,k));
                         }
                      }else{
                         ioff = 0;
                         ii = i;
                         fracx = 0.0;
                      }

                      Real tmp = (1.0-fracx)*area*bc;
                      if (cell_id_arr(i,j-1,k)>=0) {
                         mat_tmp(0, 0) = mat_tmp(0, 0) + tmp*f;
                         mat_tmp(0,-1) = mat_tmp(0,-1) - tmp*f;
                      }else if (cell_id_arr(i,j+1,k)<0 || apy_arr(i,j+1,k)==0.0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,1) = mat_tmp(0,1) + tmp*   bf2[cdir];
                      }

                      if (fracx>0.0) {
                         if (cell_id_arr(ii,j-1,k)>=0 && cell_id_arr(ii,j,k)>=0) {
                            mat_tmp(ioff,-1) = mat_tmp(ioff,-1) - fracx*area*by_arr(ii,j,k)*f;
                            mat_tmp(ioff, 0) = mat_tmp(ioff, 0) + fracx*area*by_arr(ii,j,k)*f;
                         }else if (cell_id_arr(ii,j+1,k)<0 || apy_arr(ii,j+1,k)==0.0) {
                            mat_tmp(ioff,0) = mat_tmp(ioff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,0) = mat_tmp(ioff,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,1) = mat_tmp(ioff,1) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 3;
                   f = fac[1];
                   area = apy_arr(i,j+1,k);
                   bc = by_arr(i,j+1,k);
                   if (area>0.0) {
                      if (area!=1.0) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j+1,k)));
                         ii = i+ioff;
                         if (cell_id_arr(ii,j,k)<0 && cell_id_arr(ii,j+1,k)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcy_arr(i,j+1,k));
                         }
                      }else{
                         ioff = 0;
                         ii = i;
                         fracx = 0.0;
                      }

                      Real tmp = (1.0-fracx)*area*bc;
                      if (cell_id_arr(i,j+1,k)>=0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*f;
                         mat_tmp(0,1) = mat_tmp(0,1) - tmp*f;
                      }else if (cell_id_arr(i,j-1,k)<0 || apy_arr(i,j,k)==0.0) {
                         mat_tmp(0,0) = mat_tmp(0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0, 0) = mat_tmp(0, 0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,-1) = mat_tmp(0,-1) + tmp*   bf2[cdir];
                      }

                      if (fracx>0.0) {
                         if (cell_id_arr(ii,j,k)>=0 && cell_id_arr(ii,j+1,k)>=0) {
                            mat_tmp(ioff,0) = mat_tmp(ioff,0) + fracx*area*by_arr(ii,j+1,k)*f;
                            mat_tmp(ioff,1) = mat_tmp(ioff,1) - fracx*area*by_arr(ii,j+1,k)*f;
                         }else if (cell_id_arr(ii,j-1,k)<0 || apy_arr(ii,j,k)==0.0) {
                            mat_tmp(ioff,0) = mat_tmp(ioff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff, 0) = mat_tmp(ioff, 0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,-1) = mat_tmp(ioff,-1) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   if (is_dirichlet){
                     Real anorm = std::sqrt((apx_arr(i,j,k) - apx_arr(i+1,j,k))*(apx_arr(i,j,k) - apx_arr(i+1,j,k)) 
                              + (apy_arr(i,j,k) - apy_arr(i,j+1,k))*(apy_arr(i,j,k) - apy_arr(i,j+1,k))); 

                     Real anorminv = 1.0/anorm;
                     Real anrmx = (apx_arr(i,j,k) - apx_arr(i+1,j,k))*anorminv;
                     Real anrmy = (apy_arr(i,j,k) - apy_arr(i,j+1,k))*anorminv;
                     Real sx   = Math::copysign(1.0,anrmx);
                     Real sy   = Math::copysign(1.0,anrmy);
                     Real bctx = bcen_arr(i,j,k,0);
                     Real bcty = bcen_arr(i,j,k,1);
                     Real dg   = get_dx_eb(vfrc_arr(i,j,k)) / max(Math::abs(anrmx), Math::abs(anrmy));
                     Real gx   = bctx - dg*anrmx;
                     Real gy   = bcty - dg*anrmy;
                     int ioff = -static_cast<int>(sx);
                     int joff = -static_cast<int>(sy);
                     phig1(0) = 1.0 + gx*sx + gy*sy + gx*gy*sx*sy;
                     phig1(1) =     - gx*sx         - gx*gy*sx*sy;
                     phig1(2) =             - gy*sy - gx*gy*sx*sy;
                     phig1(3) =                     + gx*gy*sx*sy;

                     for (int ii=0; ii<4; ii++){
                         feb(ii) = -phig1(ii) * (ba_arr(i,j,k) * beb_arr(i,j,k) / dg);
                     }
                     mat_tmp(0   , 0  ) = mat_tmp(0   , 0  ) - feb(0)*fac[0];
                     mat_tmp(ioff, 0  ) = mat_tmp(ioff, 0  ) - feb(1)*fac[0];
                     mat_tmp(0   ,joff) = mat_tmp(0   ,joff) - feb(2)*fac[0];
                     mat_tmp(ioff,joff) = mat_tmp(ioff,joff) - feb(3)*fac[0];

                   }

                   for (int jj=-1; jj<=1; jj++){
                       for (int ii=-1; ii<=1; ii++){
                           mat_tmp(ii,jj) = mat_tmp(ii,jj) * (1.0/vfrc_arr(i,j,k));
                       }
                   }
                   mat_tmp(0,0) = mat_tmp(0,0) + sa*a_arr(i,j,k);
            }
                   diag_arr(i,j,k) = 1.0/mat_tmp(0,0);
                   for (int jj=-1; jj<=1; jj++){
                       for (int ii=-1; ii<=1; ii++){
                           if (mat_tmp(ii,jj) != 0.0 && cell_id_arr(i+ii,j+jj,k)>=0){
                               int imat = iterPtr[1];
                               ncols[irow]  = ncols[irow] +1;
                               colsg[imat]  = cell_id_arr(i+ii,j+jj,k);
                               matg[imat] = mat_tmp(ii,jj)*diag_arr(i,j,k);
                               Gpu::Atomic::Add(&iterPtr[1], 1);
                           }
                       }
                   }
                   Gpu::Atomic::Add(&iterPtr[0], 1);
        }
}

#endif

}
#endif
