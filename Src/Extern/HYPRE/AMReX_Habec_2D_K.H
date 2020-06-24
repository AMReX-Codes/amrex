#ifndef AMREX_HABEC_2D_H_
#define AMREX_HABEC_2D_H_

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Hypre.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBABecLap_K.H>
#endif

namespace amrex {

void amrex_hpacoef (Box const& box, 
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& a,
                    Real& sa);

void amrex_hpbcoef (Box const& box,
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& b,
                    Real& sb,
                    const Real* dx,
                    int& idim);

void amrex_hpmat (Box const& box,
                  BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                  const FArrayBox& b,
                  const Mask& msk,
                  Real& sb,
                  const Real* dx,
                  int& cdir,
                  const int& bct,
                  const Real& bcl,
                  const int& bho);

void amrex_hpdiag (Box const& box,
                   BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                   FArrayBox& diag);

void amrex_hpijmatrix (Box const& box,
                       const HYPRE_Int& nrows, 
                       HYPRE_Int* ncols,
                       HYPRE_Int* rows, 
                       HYPRE_Int* colsg,
                       Real* matg, 
                       BaseFab<int>& cell_id,
                       HYPRE_Int& offset, FArrayBox& diaginv,
                       FArrayBox& acoefs, FArrayBox& bcoefsx,
                       FArrayBox& bcoefsy,
                       Real& sa, Real& sb, const Real* dx, 
                       GpuArray<int,AMREX_SPACEDIM*2>& bct,
                       GpuArray<Real,AMREX_SPACEDIM*2> bcl,
                       const int& bho);

#ifdef AMREX_USE_EB

void amrex_hpeb_fill_cellid (Box const& box,
                             int& nrows,
                             BaseFab<HYPRE_Int>& cell_id, 
                             const EBCellFlagFab& flag);

void amrex_hpeb_copy_from_vec (Box const& box,
                               FArrayBox& a,
                               Real* v,
                               const EBCellFlagFab& flag);

void amrex_hpeb_copy_to_vec (Box const& box,
                             FArrayBox& a,
                             Real* v,
                             const EBCellFlagFab& flag);

void amrex_hpeb_ijmatrix (Box const& box,
                          const HYPRE_Int& nrows,
                          HYPRE_Int* ncols,
                          HYPRE_Int* rows,
                          HYPRE_Int* colsg,
                          Real* matg,
                          BaseFab<int>& cell_id,
                          HYPRE_Int& offset, FArrayBox& diaginv,
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
                          const int& bho);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hpeb_ij (int i, int j, int k,
              const HYPRE_Int& nrows,
              HYPRE_Int* ncols,
              HYPRE_Int* rows,
              HYPRE_Int* colsg,
              Real* matg,
              const GpuArray<Real,AMREX_SPACEDIM>& fac,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bf1,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bf2,
              const GpuArray<Real,AMREX_SPACEDIM*2>& bflo,
              const Dim3 lo,
              const Dim3 hi,
              const Real& sa, const int& bho, const int& is_dirichlet,
              int* iterPtr,
              Array4<int> const& cell_id_arr,
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
        Array1D<Real,0,4> phig1;
        Array1D<Real,0,4> feb;

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
