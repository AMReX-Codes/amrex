#ifndef AMREX_HABEC_2D_H_
#define AMREX_HABEC_2D_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBABecLap_K.H>
#endif

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_mat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, int i, int j, int k,
                Dim3 const& boxlo, Dim3 const& boxhi,
                Real sa, Array4<Real const> const& a,
                Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int bho,
                GpuArray<Array4<int const>, AMREX_SPACEDIM*2> const& msk,
                Array4<Real> const& diaginv)
{
    sten[1] = -(sb / (dx[0]*dx[0])) * b[0](i,j,k);
    sten[2] = -(sb / (dx[0]*dx[0])) * b[0](i+1,j,k);
    sten[3] = -(sb / (dx[1]*dx[1])) * b[1](i,j,k);
    sten[4] = -(sb / (dx[1]*dx[1])) * b[1](i,j+1,k);
    sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4]);
    if (sa != Real(0.0)) {
        sten[0] += sa * a(i,j,k);
    }

    // xlo
    if (i == boxlo.x) {
        int cdir = Orientation(Direction::x, Orientation::low);
        if (msk[cdir](i-1,j,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
        }
    }

    // xhi
    if (i == boxhi.x) {
        int cdir = Orientation(Direction::x, Orientation::high);
        if (msk[cdir](i+1,j,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
        }
    }

    // ylo
    if (j == boxlo.y) {
        int cdir = Orientation(Direction::y, Orientation::low);
        if (msk[cdir](i,j-1,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[1](i,j,k);
            sten[3] = Real(0.0);
            sten[4] += bf2 * b[1](i,j,k);
        }
    }

    // yhi
    if (j == boxhi.y) {
        int cdir = Orientation(Direction::y, Orientation::high);
        if (msk[cdir](i,j+1,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[1](i,j+1,k);
            sten[3] += bf2 * b[1](i,j+1,k);
            sten[4] = Real(0.0);
        }
    }

    diaginv(i,j,k) = Real(1.0) / sten[0];
    sten[0] = Real(1.0);
    for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
        sten[m] *= diaginv(i,j,k);
    }
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_ijmat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, Array4<Int> const& ncols,
                  Array4<Real> const& diaginv, int i, int j, int k,
                  Array4<Int const> const& cell_id,
                  Real sa, Array4<Real const> const& a,
                  Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                  GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                  GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int bho,
                  Array4<int const> const& osm)
{
    if (!osm || osm(i,j,k) != 0) {
        sten[1] = -(sb / (dx[0]*dx[0])) * b[0](i,j,k);
        sten[2] = -(sb / (dx[0]*dx[0])) * b[0](i+1,j,k);
        sten[3] = -(sb / (dx[1]*dx[1])) * b[1](i,j,k);
        sten[4] = -(sb / (dx[1]*dx[1])) * b[1](i,j+1,k);
        sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4]);
        if (sa != Real(0.0)) {
            sten[0] += sa * a(i,j,k);
        }

        // xlo
        if (cell_id(i-1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
        }

        // xhi
        if (cell_id(i+1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
        }

        // ylo
        if (cell_id(i,j-1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[1](i,j,k);
            sten[3] = Real(0.0);
            sten[4] += bf2 * b[1](i,j,k);
        }

        // yhi
        if (cell_id(i,j+1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[1](i,j+1,k);
            sten[3] += bf2 * b[1](i,j+1,k);
            sten[4] = Real(0.0);
        }
    } else {
        sten[0] = Real(1.0);
        for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
            sten[m] = Real(0.0);
        }
    }

    diaginv(i,j,k) = Real(1.0) / sten[0];
    sten[0] = Real(1.0);
    for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
        sten[m] *= diaginv(i,j,k);
    }

    ncols(i,j,k) = 0;
    for (int m = 0; m < 2*AMREX_SPACEDIM+1; ++m) {
        ncols(i,j,k) += (sten[m] != Real(0.0));
    }
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_cols (GpuArray<Int,2*AMREX_SPACEDIM+1>& sten, int i, int j, int /*k*/,
                 Array4<Int const> const& cell_id)
{
    sten[0] = cell_id(i  ,j  ,0);
    sten[1] = cell_id(i-1,j  ,0);
    sten[2] = cell_id(i+1,j  ,0);
    sten[3] = cell_id(i  ,j-1,0);
    sten[4] = cell_id(i  ,j+1,0);
}

#ifdef AMREX_USE_EB

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_cols_eb (GpuArray<Int,9>& sten, int i, int j, int /*k*/,
                    Array4<Int const> const& cell_id)
{
    sten[0] = cell_id(i-1,j-1,0);
    sten[1] = cell_id(i  ,j-1,0);
    sten[2] = cell_id(i+1,j-1,0);
    sten[3] = cell_id(i-1,j  ,0);
    sten[4] = cell_id(i  ,j  ,0);
    sten[5] = cell_id(i+1,j  ,0);
    sten[6] = cell_id(i-1,j+1,0);
    sten[7] = cell_id(i  ,j+1,0);
    sten[8] = cell_id(i+1,j+1,0);
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_ijmat_eb (GpuArray<Real,9>& sten, Array4<Int> const& ncols,
                     Array4<Real> const& diaginv, int i, int j, int k,
                     Array4<Int const> const& cell_id,
                     Real sa, Array4<Real const> const& a,
                     Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                     GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                     GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                     GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int bho,
                     Array4<const EBCellFlag> const& flag,
                     Array4<Real const> const& vfrc,
                     Array4<Real const> const& apx,
                     Array4<Real const> const& apy,
                     Array4<Real const> const& fcx,
                     Array4<Real const> const& fcy,
                     Array4<Real const> const& ba,
                     Array4<Real const> const& bcen,
                     Array4<Real const> const& beb)
{
    for (int m = 0; m < 9; ++m) {
        sten[m] = Real(0.0);
    }

    auto& mat_tmp = reinterpret_cast<Array2D<Real,-1,1,-1,1>&>(sten);

    GpuArray<Real,AMREX_SPACEDIM> fac{sb/(dx[0]*dx[0]), sb/(dx[1]*dx[1])};

    if (flag(i,j,k).isRegular())
    {
        mat_tmp(0,0) = fac[0]*(b[0](i,j,k)+b[0](i+1,j,k))
                     + fac[1]*(b[1](i,j,k)+b[1](i,j+1,k));
        mat_tmp(-1, 0) = -fac[0]*b[0](i,j,k);
        mat_tmp( 1, 0) = -fac[0]*b[0](i+1,j,k);
        mat_tmp( 0,-1) = -fac[1]*b[1](i,j,k);
        mat_tmp( 0, 1) = -fac[1]*b[1](i,j+1,k);

        if (cell_id(i-1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0) += bf1*b[0](i,j,k);
            mat_tmp(-1,0) = Real(0.0);
            mat_tmp(1,0) += bf2*b[0](i,j,k);
        }

        if (cell_id(i+1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0) += bf1*b[0](i+1,j,k);
            mat_tmp(1,0) = Real(0.0);
            mat_tmp(-1,0) += bf2*b[0](i+1,j,k);
        }

        if (cell_id(i,j-1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0) += bf1*b[1](i,j,k);
            mat_tmp(0,-1) = Real(0.0);
            mat_tmp(0,1) += bf2*b[1](i,j,k);
        }

        if (cell_id(i,j+1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0) += bf1*b[1](i,j+1,k);
            mat_tmp(0,1) = Real(0.0);
            mat_tmp(0,-1) += bf2*b[1](i,j+1,k);
        }

        if (sa != Real(0.0)) {
            mat_tmp(0,0) += sa*a(i,j,k);
        }
    }
    else if (flag(i,j,k).isSingleValued())
    {
        int cdir = Orientation(Direction::x, Orientation::low);
        Real area = apx(i,j,k);

        if (area > Real(0.0))
        {
            int joff, jj;
            Real fracy;
            Real f = fac[0];

            if (area != Real(1.0)) {
                joff = static_cast<int>(std::copysign(Real(1.0), fcx(i,j,k)));
                jj = j+joff;
                if (cell_id(i-1,jj,k) < 0 && cell_id(i,jj,k) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcx(i,j,k));
                }
            } else {
                joff = 0;
                jj = j;
                fracy = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[0], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracy)*area*b[0](i,j,k);
            if (cell_id(i-1,j,k) >= 0) {
                mat_tmp( 0,0) += tmp*f;
                mat_tmp(-1,0) -= tmp*f;
            } else if (cell_id(i+1,j,k) < 0 || apx(i+1,j,k) == Real(0.0)) {
                mat_tmp(0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0) += tmp*(f+bf1);
                mat_tmp(1,0) += tmp*   bf2;
            }

            if (fracy > Real(0.0)) {
                if (cell_id(i-1,jj,k) >= 0 && cell_id(i,jj,k) >= 0) {
                    mat_tmp(-1,joff) -= fracy*area*b[0](i,jj,k)*f;
                    mat_tmp( 0,joff) += fracy*area*b[0](i,jj,k)*f;
                } else if (cell_id(i+1,jj,k) < 0 || apx(i+1,jj,k) == Real(0.0)) {
                    mat_tmp( 0,joff) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,joff) += tmp*(f+bf1);
                    mat_tmp( 1,joff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::x, Orientation::high);
        area = apx(i+1,j,k);
        if (area > Real(0.0))
        {
            int joff, jj;
            Real fracy;
            Real f = fac[0];

            if (area != Real(1.0)) {
                joff = static_cast<int>(std::copysign(Real(1.0), fcx(i+1,j,k)));
                jj = j+joff;
                if (cell_id(i,jj,k) < 0 && cell_id(i+1,jj,k) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcx(i+1,j,k));
                }
            } else {
                joff = 0;
                jj = j;
                fracy = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[0], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracy)*area*b[0](i+1,j,k);
            if (cell_id(i+1,j,k) >= 0) {
                mat_tmp(0,0) += tmp*f;
                mat_tmp(1,0) -= tmp*f;
            } else if (cell_id(i-1,j,k) < 0 || apx(i,j,k) == Real(0.0)) {
                mat_tmp(0,0) += tmp*(f+bflo);
            } else {
                mat_tmp( 0,0) += tmp*(f+bf1);
                mat_tmp(-1,0) += tmp*   bf2;
            }

            if (fracy > Real(0.0)) {
                if (cell_id(i,jj,k) >= 0 && cell_id(i+1,jj,k) >= 0) {
                    mat_tmp(0,joff) += fracy*area*b[0](i+1,jj,k)*f;
                    mat_tmp(1,joff) -= fracy*area*b[0](i+1,jj,k)*f;
                } else if (cell_id(i-1,jj,k) < 0 || apx(i,jj,k) == Real(0.0)) {
                    mat_tmp(0,joff) += tmp*(f+bflo);
                } else{
                    mat_tmp( 0,joff) += tmp*(f+bf1);
                    mat_tmp(-1,joff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::y, Orientation::low);
        area = apy(i,j,k);
        if (area > Real(0.0))
        {
            int ioff, ii;
            Real fracx ;
            Real f = fac[1];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j,k)));
                ii = i+ioff;
                if (cell_id(ii,j-1,k) < 0 && cell_id(ii,j,k) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcy(i,j,k));
                }
            } else {
                ioff = 0;
                ii = i;
                fracx = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[1], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*area*b[1](i,j,k);
            if (cell_id(i,j-1,k) >= 0) {
                mat_tmp(0, 0) += tmp*f;
                mat_tmp(0,-1) -= tmp*f;
            } else if (cell_id(i,j+1,k) < 0 || apy(i,j+1,k) == Real(0.0)) {
                mat_tmp(0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0) += tmp*(f+bf1);
                mat_tmp(0,1) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                if (cell_id(ii,j-1,k) >= 0 && cell_id(ii,j,k) >= 0) {
                    mat_tmp(ioff,-1) -= fracx*area*b[1](ii,j,k)*f;
                    mat_tmp(ioff, 0) += fracx*area*b[1](ii,j,k)*f;
                } else if (cell_id(ii,j+1,k) < 0 || apy(ii,j+1,k) == Real(0.0)) {
                    mat_tmp(ioff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,0) += tmp*(f+bf1);
                    mat_tmp(ioff,1) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::y, Orientation::high);
        area = apy(i,j+1,k);
        if (area > Real(0.0))
        {
            int ioff, ii;
            Real fracx;
            Real f = fac[1];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j+1,k)));
                ii = i+ioff;
                if (cell_id(ii,j,k) < 0 && cell_id(ii,j+1,k) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcy(i,j+1,k));
                }
            } else {
                ioff = 0;
                ii = i;
                fracx = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[1], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*area*b[1](i,j+1,k);
            if (cell_id(i,j+1,k) >= 0) {
                mat_tmp(0,0) += tmp*f;
                mat_tmp(0,1) -= tmp*f;
            } else if (cell_id(i,j-1,k) < 0 || apy(i,j,k) == Real(0.0)) {
                mat_tmp(0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0, 0) += tmp*(f+bf1);
                mat_tmp(0,-1) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                if (cell_id(ii,j,k) >= 0 && cell_id(ii,j+1,k) >= 0) {
                    mat_tmp(ioff,0) += fracx*area*b[1](ii,j+1,k)*f;
                    mat_tmp(ioff,1) -= fracx*area*b[1](ii,j+1,k)*f;
                } else if (cell_id(ii,j-1,k) < 0 || apy(ii,j,k) == Real(0.0)) {
                    mat_tmp(ioff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff, 0) += tmp*(f+bf1);
                    mat_tmp(ioff,-1) += tmp*   bf2;
                }
            }
        }

        if (beb) {
            Real anorm = std::sqrt((apx(i,j,k) - apx(i+1,j,k))*(apx(i,j,k) - apx(i+1,j,k))
                                 + (apy(i,j,k) - apy(i,j+1,k))*(apy(i,j,k) - apy(i,j+1,k)));

            Real anorminv = Real(1.0)/anorm;
            Real anrmx = (apx(i,j,k) - apx(i+1,j,k))*anorminv;
            Real anrmy = (apy(i,j,k) - apy(i,j+1,k))*anorminv;
            Real sx   = std::copysign(Real(1.0),anrmx);
            Real sy   = std::copysign(Real(1.0),anrmy);
            Real bctx = bcen(i,j,k,0);
            Real bcty = bcen(i,j,k,1);
            Real dg   = get_dx_eb(vfrc(i,j,k)) / amrex::max(std::abs(anrmx), std::abs(anrmy));
            Real gx   = bctx - dg*anrmx;
            Real gy   = bcty - dg*anrmy;
            int ioff = -static_cast<int>(sx);
            int joff = -static_cast<int>(sy);
            Array1D<Real,0,3> phig1{Real(1.0) + gx*sx + gy*sy + gx*gy*sx*sy,
                                              - gx*sx         - gx*gy*sx*sy,
                                                      - gy*sy - gx*gy*sx*sy,
                                                              + gx*gy*sx*sy};

            Array1D<Real,0,3> feb;
            for (int ii=0; ii<4; ii++){
                feb(ii) = -phig1(ii) * (ba(i,j,k) * beb(i,j,k) / dg);
            }
            mat_tmp(0   , 0  ) -= feb(0)*fac[0];
            mat_tmp(ioff, 0  ) -= feb(1)*fac[0];
            mat_tmp(0   ,joff) -= feb(2)*fac[0];
            mat_tmp(ioff,joff) -= feb(3)*fac[0];
        }

        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            mat_tmp(ii,jj) *= Real(1.0)/vfrc(i,j,k);
        }}

        if (sa != Real(0.0)) {
            mat_tmp(0,0) += sa*a(i,j,k);
        }
    }
    else
    {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            mat_tmp(ii,jj) = Real(0.0);
        }}
        mat_tmp(0,0) = Real(1.0);
    }

    diaginv(i,j,k) = Real(1.0) / mat_tmp(0,0);
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        mat_tmp(ii,jj) *= diaginv(i,j,k);
    }}
    mat_tmp(0,0) = Real(1.0);

    ncols(i,j,k) = 0;
    for (int m = 0; m < 9; ++m) {
        ncols(i,j,k) += (sten[m] != Real(0.0));
    }
}

#endif

}
#endif
