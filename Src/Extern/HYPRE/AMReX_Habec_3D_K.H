#ifndef AMREX_HABEC_3D_H_
#define AMREX_HABEC_3D_H_

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
                       FArrayBox& bcoefsy, FArrayBox& bcoefsz,
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
                          FArrayBox& bcoefsy, FArrayBox& bcoefsz,
                          const EBCellFlagFab& flag,
                          const FArrayBox& vfrc,
                          const FArrayBox& apx, const FArrayBox& apy,
                          const FArrayBox& apz, const FArrayBox& fcx,
                          const FArrayBox& fcy, const FArrayBox& fcz,
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
              Array4<Real const> const& bz_arr,
              Array4<Real> const& diag_arr,
              Array4<const EBCellFlag> const& flag_arr,
              Array4<Real const> const& apx_arr,
              Array4<Real const> const& fcx_arr,
              Array4<Real const> const& apy_arr,
              Array4<Real const> const& fcy_arr,
              Array4<Real const> const& apz_arr,
              Array4<Real const> const& fcz_arr,
              Array4<Real const> const& bcen_arr,
              Array4<Real const> const& vfrc_arr,
              Array4<Real const> const& ba_arr,
              Array4<Real const> const& beb_arr) noexcept
{
        Array3D<Real,-1,1,-1,1,-1,1> mat_tmp;
        Array1D<Real,0,7> phig1;
        Array1D<Real,0,7> feb;

        if (flag_arr(i,j,k).isCovered()){
            diag_arr(i,j,k) = 0.0;
        }else{
            int irow = iterPtr[0];
            rows[irow]  = cell_id_arr(i,j,k);
            ncols[irow] = 0;
            for (int kk=-1; kk<=1; kk++){
                for (int jj=-1; jj<=1; jj++){
                    for (int ii=-1; ii<=1; ii++){
                         mat_tmp(ii,jj,kk) = 0.0;
                    }
                }
            }

            if (flag_arr(i,j,k).isRegular()){
                mat_tmp(0,0,0) = sa*a_arr(i,j,k) + fac[0]*(bx_arr(i,j,k)+bx_arr(i+1,j,k))
                                                 + fac[1]*(by_arr(i,j,k)+by_arr(i,j+1,k))
                                                 + fac[2]*(bz_arr(i,j,k)+bz_arr(i,j,k+1));
                mat_tmp(-1, 0, 0) = -fac[0]*bx_arr(i,j,k);
                mat_tmp( 1, 0, 0) = -fac[0]*bx_arr(i+1,j,k);
                mat_tmp( 0,-1, 0) = -fac[1]*by_arr(i,j,k);
                mat_tmp( 0, 1, 0) = -fac[1]*by_arr(i,j+1,k);
                mat_tmp( 0, 0,-1) = -fac[2]*bz_arr(i,j,k);
                mat_tmp( 0, 0, 1) = -fac[2]*bz_arr(i,j,k+1);

                if (i == lo.x && cell_id_arr(i-1,j,k)<0) {
                   mat_tmp(0,0,0)  = mat_tmp(0,0,0) + bf1[0]*bx_arr(i,j,k);
                   mat_tmp(-1,0,0) = 0.0;
                   mat_tmp(1,0,0)  = mat_tmp(1,0,0) + bf2[0]*bx_arr(i,j,k);
                }

                if (i == hi.x && cell_id_arr(i+1,j,k)<0) {
                   mat_tmp(0,0,0)  = mat_tmp(0,0,0) + bf1[3]*bx_arr(i+1,j,k);
                   mat_tmp(1,0,0)  = 0.0;
                   mat_tmp(-1,0,0) = mat_tmp(-1,0,0) + bf2[3]*bx_arr(i+1,j,k);
                }

                if (j == lo.y && cell_id_arr(i,j-1,k)<0) {
                   mat_tmp(0,0,0)  = mat_tmp(0,0,0) + bf1[1]*by_arr(i,j,k);
                   mat_tmp(0,-1,0) = 0.0;
                   mat_tmp(0,1,0)  = mat_tmp(0,1,0) + bf2[1]*by_arr(i,j,k);
                }

                if (j == hi.y && cell_id_arr(i,j+1,k)<0) {
                   mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1[4]*by_arr(i,j+1,k);
                   mat_tmp(0,1,0) = 0.0;
                   mat_tmp(0,-1,0) = mat_tmp(0,-1,0) + bf2[4]*by_arr(i,j+1,k);
                }

                if (k == lo.z && cell_id_arr(i,j,k-1)<0) {
                   mat_tmp(0,0,0)  = mat_tmp(0,0,0) + bf1[2]*bz_arr(i,j,k);
                   mat_tmp(0,0,-1) = 0.0;
                   mat_tmp(0,0,1)  = mat_tmp(0,0,1) + bf2[2]*bz_arr(i,j,k);
                }

                if (k == hi.z && cell_id_arr(i,j,k+1)<0) {
                   mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1[5]*bz_arr(i,j,k+1);
                   mat_tmp(0,0,1) = 0.0;
                   mat_tmp(0,0,-1) = mat_tmp(0,0,-1) + bf2[5]*bz_arr(i,j,k+1);
                }
            }else{
                   int cdir = 0;
                   int joff = 0;
                   int koff = 0;
                   int jj = 0;
                   int kk = 0;
                   Real f = fac[0];
                   Real area = apx_arr(i,j,k);
                   Real abc = area*bx_arr(i,j,k);
                   Real fracy = 0.0;
                   Real fracz = 0.0;
                   if (area > 0.0) {
                      if (area != 1.0) {
                         joff = static_cast<int>(Math::copysign(1.0, fcx_arr(i,j,k,0)));
                         koff = static_cast<int>(Math::copysign(1.0, fcx_arr(i,j,k,1)));  
                         jj = j+joff;
                         kk = k+koff;
                         if (cell_id_arr(i-1,jj,k)<0 && cell_id_arr(i,jj,k)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcx_arr(i,j,k,0));
                         }
                         if (cell_id_arr(i-1,j,kk)<0 && cell_id_arr(i,j,kk)<0) {
                            fracz = 0.0;
                         }else{
                            fracz = Math::abs(fcx_arr(i,j,k,1));
                         }
                         if (cell_id_arr(i-1,jj,kk)<0 && cell_id_arr(i,jj,kk)<0 && (fracy*fracz)>0.0) {
                            fracy = 0.0;
                            fracz = 0.0;
                         }
                      }else{
                         joff = 0;
                         koff = 0;
                         jj = j;
                         kk = k;
                         fracy = 0.0;
                         fracz = 0.0;
                      }

                      Real tmp = (1.0-fracy)*(1.0-fracz)*abc;
                      // cell(i-1,j,k) is not covered because area > 0
                      if (cell_id_arr(i-1,j,k)>=0) {
                         mat_tmp( 0,0,0) = mat_tmp( 0,0,0) + tmp*f;
                         mat_tmp(-1,0,0) = mat_tmp(-1,0,0) - tmp*f;
                      }else if (cell_id_arr(i+1,j,k)<0 || apx_arr(i+1,j,k)==0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(1,0,0) = mat_tmp(1,0,0) + tmp*   bf2[cdir];
                      }

                      if (fracy>0.0) {
                         tmp = fracy*(1.0-fracz)*area*bx_arr(i,jj,k);
                         if (cell_id_arr(i-1,jj,k)>=0 && cell_id_arr(i,jj,k)>=0) {
                            mat_tmp(-1,joff,0) = mat_tmp(-1,joff,0) - tmp*f;
                            mat_tmp( 0,joff,0) = mat_tmp( 0,joff,0) + tmp*f;
                         }else if (cell_id_arr(i+1,jj,k)<0 || apx_arr(i+1,jj,k)==0 || bho==0) {
                            mat_tmp( 0,joff,0) = mat_tmp( 0,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff,0) = mat_tmp( 0,joff,0) + tmp*(f+bf1[cdir]);
                            mat_tmp( 1,joff,0) = mat_tmp( 1,joff,0) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracz>0.0) {
                         tmp = fracz*(1.0-fracy)*area*bx_arr(i,j,kk);
                         if (cell_id_arr(i-1,j,kk)>=0 && cell_id_arr(i,j,kk)>=0) {
                            mat_tmp(-1,0,koff) = mat_tmp(-1,0,koff) - tmp*f;
                            mat_tmp( 0,0,koff) = mat_tmp( 0,0,koff) + tmp*f;
                         }else if (cell_id_arr(i+1,j,kk)<0 || apx_arr(i+1,j,kk)==0 || bho==0) {
                            mat_tmp( 0,0,koff) = mat_tmp( 0,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,0,koff) = mat_tmp( 0,0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp( 1,0,koff) = mat_tmp( 1,0,koff) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracy>0.0 && fracz>0.0) {
                         tmp = fracy*fracz*area*bx_arr(i,jj,kk);
                         if (cell_id_arr(i-1,jj,kk)>=0 && cell_id_arr(i,jj,kk)>=0) {
                            mat_tmp(-1,joff,koff) = mat_tmp(-1,joff,koff) - tmp*f;
                            mat_tmp( 0,joff,koff) = mat_tmp( 0,joff,koff) + tmp*f;
                         }else if (cell_id_arr(i+1,jj,kk)<0 || apx_arr(i+1,jj,kk)==0 || bho==0) {
                            mat_tmp( 0,joff,koff) = mat_tmp( 0,joff,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff,koff) = mat_tmp( 0,joff,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp( 1,joff,koff) = mat_tmp( 1,joff,koff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 3;
                   f = fac[0];
                   area = apx_arr(i+1,j,k);
                   abc = area*bx_arr(i+1,j,k);
                   if (area>0.0) {
                      if (area!=1.0) {
                         joff = static_cast<int>(Math::copysign(1.0, fcx_arr(i+1,j,k,0)));
                         koff = static_cast<int>(Math::copysign(1.0, fcx_arr(i+1,j,k,1)));
                         jj = j+joff;
                         kk = k+koff;
                         if (cell_id_arr(i,jj,k)<0 && cell_id_arr(i+1,jj,k)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcx_arr(i+1,j,k,0));
                         }
                         if (cell_id_arr(i,j,kk)<0 && cell_id_arr(i+1,j,kk)<0) {
                            fracz = 0.0;
                         }else{
                            fracz = Math::abs(fcx_arr(i+1,j,k,1));
                         }
                         if (cell_id_arr(i,jj,kk)<0 && cell_id_arr(i+1,jj,kk)<0 && (fracy*fracz)>0.0) {
                            fracy = 0.0;
                            fracz = 0.0;
                         }
                      }else{
                         joff = 0;
                         koff = 0;
                         jj = j;
                         kk = k;
                         fracy = 0.0;
                         fracz = 0.0;
                      }

                      Real tmp = (1.0-fracy)*(1.0-fracz)*abc;
                      if (cell_id_arr(i+1,j,k)>=0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*f;
                         mat_tmp(1,0,0) = mat_tmp(1,0,0) - tmp*f;
                      }else if (cell_id_arr(i-1,j,k)<0 || apx_arr(i-1,j,k)==0.0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp( 0,0,0) = mat_tmp( 0,0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(-1,0,0) = mat_tmp(-1,0,0) + tmp*   bf2[cdir];
                      }

                      if (fracy>0.0) {
                         tmp = fracy*(1.0-fracz)*area*bx_arr(i+1,jj,k);
                         if (cell_id_arr(i,jj,k)>=0 && cell_id_arr(i+1,jj,k)>=0) {
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*f;
                            mat_tmp(1,joff,0) = mat_tmp(1,joff,0) - tmp*f;
                         }else if (cell_id_arr(i-1,jj,k)<0 || apx_arr(i-1,jj,k)==0.0 || bho==0) {
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff,0) = mat_tmp( 0,joff,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(-1,joff,0) = mat_tmp(-1,joff,0) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracz>0.0) {
                         tmp = fracz*(1.0-fracy)*area*bx_arr(i+1,j,kk);
                         if (cell_id_arr(i,j,kk)>=0 && cell_id_arr(i+1,j,kk)>=0) {
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*f;
                            mat_tmp(1,0,koff) = mat_tmp(1,0,koff) - tmp*f;
                         }else if (cell_id_arr(i-1,j,kk)<0 || apx_arr(i-1,j,kk)==0 || bho==0) {
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,0,koff) = mat_tmp( 0,0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(-1,0,koff) = mat_tmp(-1,0,koff) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracy>0.0 && fracz>0.0) {
                         tmp = fracy*fracz*area*bx_arr(i+1,jj,kk);
                         if (cell_id_arr(i,jj,kk)>=0 && cell_id_arr(i+1,jj,kk)>=0) {
                            mat_tmp(0,joff,koff) = mat_tmp(0,joff,koff) + tmp*f;
                            mat_tmp(1,joff,koff) = mat_tmp(1,joff,koff) - tmp*f;
                         }else if (cell_id_arr(i-1,jj,kk)<0 || apx_arr(i-1,jj,kk)==0 || bho==0) {
                            mat_tmp(0,joff,koff) = mat_tmp(0,joff,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp( 0,joff,koff) = mat_tmp( 0,joff,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(-1,joff,koff) = mat_tmp(-1,joff,koff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 1;
                   Real fracx = 0.0;
                   int ii = 0;
                   int ioff = 0;

                   f = fac[1];
                   area = apy_arr(i,j,k);
                   abc = area*by_arr(i,j,k);
                   if (area>0) {
                      if (area!=1) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j,k,0)));
                         koff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j,k,1)));
                         ii = i+ioff;
                         kk = k+koff;
                         if (cell_id_arr(ii,j-1,k)<0 && cell_id_arr(ii,j,k)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcy_arr(i,j,k,0));
                         }
                         if (cell_id_arr(i,j-1,kk)<0 && cell_id_arr(i,j,kk)<0) {
                            fracz = 0.0;
                         }else{
                            fracz = Math::abs(fcy_arr(i,j,k,1));
                         }
                         if (cell_id_arr(ii,j-1,kk)<0 && cell_id_arr(ii,j,kk)<0 && fracx*fracz>0.0) {
                            fracx = 0.0;
                            fracz = 0.0;
                         }
                      }else{
                         ioff = 0;
                         koff = 0;
                         ii = i;
                         kk = k;
                         fracx = 0.0;
                         fracz = 0.0;
                      }

                      Real tmp = (1.0-fracx)*(1.0-fracz)*abc;
                      if (cell_id_arr(i,j-1,k)>=0) {
                         mat_tmp(0, 0,0) = mat_tmp(0, 0,0) + tmp*f;
                         mat_tmp(0,-1,0) = mat_tmp(0,-1,0) - tmp*f;
                      }else if (cell_id_arr(i,j+1,k)<0 || apy_arr(i,j+1,k)==0.0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,1,0) = mat_tmp(0,1,0) + tmp*   bf2[cdir];
                      }

                      if (fracx>0.0) {
                         tmp = fracx*(1.0-fracz)*area*by_arr(ii,j,k);
                         if (cell_id_arr(ii,j-1,k)>=0 && cell_id_arr(ii,j,k)>=0) {
                            mat_tmp(ioff,-1,0) = mat_tmp(ioff,-1,0) - tmp*f;
                            mat_tmp(ioff, 0,0) = mat_tmp(ioff, 0,0) + tmp*f;
                         }else if (cell_id_arr(ii,j+1,k)<0 || apy_arr(ii,j+1,k)==0.0 || bho==0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,1,0) = mat_tmp(ioff,1,0) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracz>0.0) {
                         tmp = fracz*(1.0-fracx)*area*by_arr(i,j,kk);
                         if (cell_id_arr(i,j-1,kk)>=0 && cell_id_arr(i,j,kk)>=0) {
                            mat_tmp(0,-1,koff) = mat_tmp(0,-1,koff) - tmp*f;
                            mat_tmp(0, 0,koff) = mat_tmp(0, 0,koff) + tmp*f;
                         }else if (cell_id_arr(i,j+1,kk)<0 || apy_arr(i,j+1,kk)==0.0 || bho==0) {
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(0,1,koff) = mat_tmp(0,1,koff) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracx>0.0 && fracz>0.0) {
                         tmp = fracx*fracz*area*by_arr(ii,j,kk);
                         if (cell_id_arr(ii,j-1,kk)>=0 && cell_id_arr(ii,j,kk)>=0) {
                            mat_tmp(ioff,-1,koff) = mat_tmp(ioff,-1,koff) - tmp*f;
                            mat_tmp(ioff, 0,koff) = mat_tmp(ioff, 0,koff) + tmp*f;
                         }else if (cell_id_arr(ii,j+1,kk)<0 || apy_arr(ii,j+1,kk)==0.0 || bho==0) {
                            mat_tmp(ioff,0,koff) = mat_tmp(ioff,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,0,koff) = mat_tmp(ioff,0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,1,koff) = mat_tmp(ioff,1,koff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 4;
                   f = fac[1];
                   area = apy_arr(i,j+1,k);
                   abc = area*by_arr(i,j+1,k);
                   if (area>0.0) {
                      if (area!=1.0) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j+1,k,0)));
                         koff = static_cast<int>(Math::copysign(1.0, fcy_arr(i,j+1,k,1)));
                         ii = i+ioff;
                         kk = k+koff;
                         if (cell_id_arr(ii,j,k)<0 && cell_id_arr(ii,j+1,k)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcy_arr(i,j+1,k,0));
                         }
                         if (cell_id_arr(i,j,kk)<0 && cell_id_arr(i,j+1,kk)<0) {
                            fracz = 0.0;
                         }else{
                            fracz = Math::abs(fcy_arr(i,j+1,k,1));
                         }
                         if (cell_id_arr(ii,j,kk)<0 && cell_id_arr(ii,j+1,kk)<0 && fracx*fracz>0.0) {
                            fracx = 0.0;
                            fracz = 0.0;
                         }
                      }else{
                         ioff = 0;
                         koff = 0;
                         ii = i;
                         kk = k;
                         fracx = 0.0;
                         fracz = 0.0;
                      }

                      Real tmp = (1.0-fracx)*(1.0-fracz)*abc;
                      if (cell_id_arr(i,j+1,k)>=0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*f;
                         mat_tmp(0,1,0) = mat_tmp(0,1,0) - tmp*f;
                      }else if (cell_id_arr(i,j-1,k)<0 || apy_arr(i,j-1,k)==0.0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0, 0,0) = mat_tmp(0, 0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,-1,0) = mat_tmp(0,-1,0) + tmp*   bf2[cdir];
                      }

                      if (fracx>0.0) {
                         tmp = fracx*(1.0-fracz)*area*by_arr(ii,j+1,k);
                         if (cell_id_arr(ii,j,k)>=0 && cell_id_arr(ii,j+1,k)>=0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*f;
                            mat_tmp(ioff,1,0) = mat_tmp(ioff,1,0) - tmp*f;
                         }else if (cell_id_arr(ii,j-1,k)<0 || apy_arr(ii,j-1,k)==0.0 || bho==0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff, 0,0) = mat_tmp(ioff, 0,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,-1,0) = mat_tmp(ioff,-1,0) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracz>0.0) {
                         tmp = fracz*(1.0-fracx)*area*by_arr(i,j+1,kk);
                         if (cell_id_arr(i,j,kk)>=0 && cell_id_arr(i,j+1,kk)>=0) {
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*f;
                            mat_tmp(0,1,koff) = mat_tmp(0,1,koff) - tmp*f;
                         }else if (cell_id_arr(i,j-1,kk)<0 || apy_arr(i,j-1,kk)==0.0 || bho==0) {
                            mat_tmp(0,0,koff) = mat_tmp(0,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(0, 0,koff) = mat_tmp(0, 0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(0,-1,koff) = mat_tmp(0,-1,koff) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracx>0.0 && fracz>0.0) {
                         tmp = fracx*fracz*area*by_arr(ii,j+1,kk);
                         if (cell_id_arr(ii,j,kk)>=0 && cell_id_arr(ii,j+1,kk)>=0) {
                            mat_tmp(ioff,1,koff) = mat_tmp(ioff,1,koff) - tmp*f;
                            mat_tmp(ioff,0,koff) = mat_tmp(ioff,0,koff) + tmp*f;
                         }else if (cell_id_arr(ii,j-1,kk)<0 || apy_arr(ii,j-1,kk)==0.0 || bho==0) {
                            mat_tmp(ioff,0,koff) = mat_tmp(ioff,0,koff) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff, 0,koff) = mat_tmp(ioff, 0,koff) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,-1,koff) = mat_tmp(ioff,-1,koff) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 2;
                   f = fac[2];
                   area = apz_arr(i,j,k);
                   abc = area*bz_arr(i,j,k);
                   if (area>0.0) {
                      if (area!=1.0) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcz_arr(i,j,k,0)));
                         joff = static_cast<int>(Math::copysign(1.0, fcz_arr(i,j,k,1)));
                         ii = i+ioff;
                         jj = j+joff;
                         if (cell_id_arr(ii,j,k-1)<0 && cell_id_arr(ii,j,k)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcz_arr(i,j,k,0));
                         }
                         if (cell_id_arr(i,jj,k-1)<0 && cell_id_arr(i,jj,k)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcz_arr(i,j,k,1));
                         }
                         if (cell_id_arr(ii,jj,k-1)<0 && cell_id_arr(ii,jj,k)<0 && fracx*fracy>0.0) {
                            fracx = 0.0;
                            fracy = 0.0;
                         }
                      }else{
                         ioff = 0;
                         joff = 0;
                         ii = i;
                         jj = j;
                         fracx = 0.0;
                         fracy = 0.0;
                      }
                      Real tmp = (1.0-fracx)*(1.0-fracy)*abc;
                      if (cell_id_arr(i,j,k-1)>=0) {
                         mat_tmp(0,0, 0) = mat_tmp(0,0, 0) + tmp*f;
                         mat_tmp(0,0,-1) = mat_tmp(0,0,-1) - tmp*f;
                      }else if (cell_id_arr(i,j,k+1)<0 || apz_arr(i,j,k+1)==0.0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,0,1) = mat_tmp(0,0,1) + tmp*   bf2[cdir];
                      }
                      if (fracx>0.0) {
                         tmp = fracx*(1.0-fracy)*area*bz_arr(ii,j,k);
                         if (cell_id_arr(ii,j,k-1)>=0 && cell_id_arr(ii,j,k)>=0) {
                            mat_tmp(ioff,0,-1) = mat_tmp(ioff,0,-1) - tmp*f;
                            mat_tmp(ioff,0, 0) = mat_tmp(ioff,0, 0) + tmp*f;
                         }else if (cell_id_arr(ii,j,k+1)<0 || apz_arr(ii,j,k+1)==0.0 || bho==0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,0,1) = mat_tmp(ioff,0,1) + tmp*   bf2[cdir];
                         }
                      }
                      if (fracy>0.0) {
                         tmp = fracy*(1.0-fracx)*area*bz_arr(i,jj,k);
                         if (cell_id_arr(i,jj,k-1)>=0 && cell_id_arr(i,jj,k)>=0) {
                            mat_tmp(0,joff,-1) = mat_tmp(0,joff,-1) - tmp*f;
                            mat_tmp(0,joff, 0) = mat_tmp(0,joff, 0) + tmp*f;
                         }else if (cell_id_arr(i,jj,k+1)<0 || apz_arr(i,jj,k+1)==0.0 || bho==0) {
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(0,joff,1) = mat_tmp(0,joff,1) + tmp*   bf2[cdir];
                         }
                      }
                      if (fracx>0.0 && fracy>0.0) {
                         tmp = fracx*fracy*area*bz_arr(ii,jj,k);
                         if (cell_id_arr(ii,jj,k-1)>=0 && cell_id_arr(ii,jj,k)>=0) {
                            mat_tmp(ioff,joff,-1) = mat_tmp(ioff,joff,-1) - tmp*f;
                            mat_tmp(ioff,joff, 0) = mat_tmp(ioff,joff, 0) + tmp*f;
                         }else if (cell_id_arr(ii,jj,k+1)<0 || apz_arr(ii,jj,k+1)==0.0 || bho==0) {
                            mat_tmp(ioff,joff,0) = mat_tmp(ioff,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,joff,0) = mat_tmp(ioff,joff,0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,joff,1) = mat_tmp(ioff,joff,1) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   cdir = 5;
                   f = fac[2];
                   area = apz_arr(i,j,k+1);
                   abc = area*bz_arr(i,j,k+1);
                   if (area>0.0) {
                      if (area!=1.0) {
                         ioff = static_cast<int>(Math::copysign(1.0, fcz_arr(i,j,k+1,0)));
                         joff = static_cast<int>(Math::copysign(1.0, fcz_arr(i,j,k+1,1)));
                         ii = i+ioff;
                         jj = j+joff;
                         if (cell_id_arr(ii,j,k)<0 && cell_id_arr(ii,j,k+1)<0) {
                            fracx = 0.0;
                         }else{
                            fracx = Math::abs(fcz_arr(i,j,k+1,0));
                         }
                         if (cell_id_arr(i,jj,k)<0 && cell_id_arr(i,jj,k+1)<0) {
                            fracy = 0.0;
                         }else{
                            fracy = Math::abs(fcz_arr(i,j,k+1,1));
                         }
                         if (cell_id_arr(ii,jj,k)<0 && cell_id_arr(ii,jj,k+1)<0 && fracx*fracy>0.0) {
                            fracx = 0.0;
                            fracy = 0.0;
                         }
                      }else{
                         ioff = 0;
                         joff = 0;
                         ii = i;
                         jj = j;
                         fracx = 0.0;
                         fracy = 0.0;
                      }

                      Real tmp = (1.0-fracx)*(1.0-fracy)*abc;
                      if (cell_id_arr(i,j,k+1)>=0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*f;
                         mat_tmp(0,0,1) = mat_tmp(0,0,1) - tmp*f;
                      }else if (cell_id_arr(i,j,k-1)<0 || apz_arr(i,j,k-1)==0.0 || bho==0) {
                         mat_tmp(0,0,0) = mat_tmp(0,0,0) + tmp*(f+bflo[cdir]);
                      }else{
                         mat_tmp(0,0, 0) = mat_tmp(0,0, 0) + tmp*(f+bf1[cdir]);
                         mat_tmp(0,0,-1) = mat_tmp(0,0,-1) + tmp*   bf2[cdir];
                      }

                      if (fracx>0.0) {
                         tmp = fracx*(1.0-fracy)*area*bz_arr(ii,j,k+1);
                         if (cell_id_arr(ii,j,k)>=0 && cell_id_arr(ii,j,k+1)>=0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*f;
                            mat_tmp(ioff,0,1) = mat_tmp(ioff,0,1) - tmp*f;
                         }else if (cell_id_arr(ii,j,k-1)<0 || apz_arr(ii,j,k-1)==0.0 || bho==0) {
                            mat_tmp(ioff,0,0) = mat_tmp(ioff,0,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,0, 0) = mat_tmp(ioff,0, 0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,0,-1) = mat_tmp(ioff,0,-1) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracy>0.0) {
                         tmp = fracy*(1.0-fracx)*area*bz_arr(i,jj,k+1);
                         if (cell_id_arr(i,jj,k)>=0 && cell_id_arr(i,jj,k+1)>=0) {
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*f;
                            mat_tmp(0,joff,1) = mat_tmp(0,joff,1) - tmp*f;
                         }else if (cell_id_arr(i,jj,k-1)<0 || apz_arr(i,jj,k-1)==0.0 || bho==0) {
                            mat_tmp(0,joff,0) = mat_tmp(0,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(0,joff, 0) = mat_tmp(0,joff, 0) + tmp*(f+bf1[cdir]);
                            mat_tmp(0,joff,-1) = mat_tmp(0,joff,-1) + tmp*   bf2[cdir];
                         }
                      }

                      if (fracx>0.0 && fracy>0.0) {
                         tmp = fracx*fracy*area*bz_arr(ii,jj,k+1);
                         if (cell_id_arr(ii,jj,k+1)>=0 && cell_id_arr(ii,jj,k)>=0) {
                            mat_tmp(ioff,joff,1) = mat_tmp(ioff,joff,1) - tmp*f;
                            mat_tmp(ioff,joff,0) = mat_tmp(ioff,joff,0) + tmp*f;
                         }else if (cell_id_arr(ii,jj,k-1)<0 || apz_arr(ii,jj,k-1)==0.0 || bho==0) {
                            mat_tmp(ioff,joff,0) = mat_tmp(ioff,joff,0) + tmp*(f+bflo[cdir]);
                         }else{
                            mat_tmp(ioff,joff, 0) = mat_tmp(ioff,joff, 0) + tmp*(f+bf1[cdir]);
                            mat_tmp(ioff,joff,-1) = mat_tmp(ioff,joff,-1) + tmp*   bf2[cdir];
                         }
                      }
                   }

                   if (is_dirichlet){
                     Real anorm = std::sqrt((apx_arr(i,j,k) - apx_arr(i+1,j,k))*(apx_arr(i,j,k) - apx_arr(i+1,j,k)) 
                              + (apy_arr(i,j,k) - apy_arr(i,j+1,k))*(apy_arr(i,j,k) - apy_arr(i,j+1,k)) 
                              + (apz_arr(i,j,k) - apz_arr(i,j,k+1))*(apz_arr(i,j,k) - apz_arr(i,j,k+1)));

                     Real anorminv = 1.0/anorm;
                     Real anrmx = (apx_arr(i,j,k) - apx_arr(i+1,j,k))*anorminv;
                     Real anrmy = (apy_arr(i,j,k) - apy_arr(i,j+1,k))*anorminv;
                     Real anrmz = (apz_arr(i,j,k) - apz_arr(i,j,k+1))*anorminv;
                     Real sx   = Math::copysign(1.0,anrmx);
                     Real sy   = Math::copysign(1.0,anrmy);
                     Real sz   = Math::copysign(1.0,anrmz);
                     Real bctx = bcen_arr(i,j,k,0);
                     Real bcty = bcen_arr(i,j,k,1);
                     Real bctz = bcen_arr(i,j,k,2);
                     Real dg   = get_dx_eb(vfrc_arr(i,j,k)) / max(Math::abs(anrmx), Math::abs(anrmy), Math::abs(anrmz));
                     Real gx   = bctx - dg*anrmx;
                     Real gy   = bcty - dg*anrmy;
                     Real gz   = bctz - dg*anrmz;
                     int ioff = -static_cast<int>(sx);
                     int joff = -static_cast<int>(sy);
                     int koff = -static_cast<int>(sz);
                     gx   = sx*gx;
                     gy   = sy*gy;
                     gz   = sz*gz;
                     Real gxy  = gx*gy;
                     Real gxz  = gx*gz;
                     Real gyz  = gy*gz;
                     Real gxyz = gx*gy*gz;
                     phig1(0) = 1.0 + gx + gy + gz + gxy + gxz + gyz + gxyz;
                     phig1(1) =     - gx           - gxy - gxz       - gxyz;
                     phig1(2) =          - gy      - gxy       - gyz - gxyz;
                     phig1(3) =               - gz       - gxz - gyz - gxyz;
                     phig1(4) =                    + gxy             + gxyz;
                     phig1(5) =                          + gxz       + gxyz;
                     phig1(6) =                                + gyz + gxyz;
                     phig1(7) =                                      - gxyz;
                     for (int ii=0; ii<8; ii++){
                         feb(ii) = -phig1(ii) * (ba_arr(i,j,k) * beb_arr(i,j,k) / dg);
                     }
                     mat_tmp(0   , 0  , 0  ) = mat_tmp(0   , 0 , 0   ) - feb(0)*fac[0];
                     mat_tmp(ioff, 0  , 0  ) = mat_tmp(ioff, 0 , 0   ) - feb(1)*fac[0];
                     mat_tmp(0   ,joff, 0  ) = mat_tmp(0   ,joff, 0  ) - feb(2)*fac[0];
                     mat_tmp(0   , 0  ,koff) = mat_tmp(0   , 0  ,koff) - feb(3)*fac[0];
                     mat_tmp(ioff,joff, 0  ) = mat_tmp(ioff,joff, 0  ) - feb(4)*fac[0];
                     mat_tmp(ioff, 0  ,koff) = mat_tmp(ioff, 0  ,koff) - feb(5)*fac[0];
                     mat_tmp(0   ,joff,koff) = mat_tmp(0   ,joff,koff) - feb(6)*fac[0];
                     mat_tmp(ioff,joff,koff) = mat_tmp(ioff,joff,koff) - feb(7)*fac[0];

                   }

                   for (int kk=-1; kk<=1; kk++){
                       for (int jj=-1; jj<=1; jj++){
                           for (int ii=-1; ii<=1; ii++){
                               mat_tmp(ii,jj,kk) = mat_tmp(ii,jj,kk) * (1.0/vfrc_arr(i,j,k));
                           }
                       }
                   }
                   mat_tmp(0,0,0) = mat_tmp(0,0,0) + sa*a_arr(i,j,k);
            }
                   diag_arr(i,j,k) = 1.0/mat_tmp(0,0,0);
                   for (int kk=-1; kk<=1; kk++){
                       for (int jj=-1; jj<=1; jj++){
                           for (int ii=-1; ii<=1; ii++){
                               if (mat_tmp(ii,jj,kk) != 0.0 && cell_id_arr(i+ii,j+jj,k+kk)>=0){
                                   int imat = iterPtr[1];
                                   ncols[irow]  = ncols[irow] +1;
                                   colsg[imat]  = cell_id_arr(i+ii,j+jj,k+kk);
                                   matg[imat] = mat_tmp(ii,jj,kk)*diag_arr(i,j,k);
                                   Gpu::Atomic::Add(&iterPtr[1], 1);
                               }
                           }
                       }
                   }
                   Gpu::Atomic::Add(&iterPtr[0], 1);
        }
}

#endif

}
#endif
