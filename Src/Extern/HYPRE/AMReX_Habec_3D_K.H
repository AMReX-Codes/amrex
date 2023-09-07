#ifndef AMREX_HABEC_3D_H_
#define AMREX_HABEC_3D_H_
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
    sten[5] = -(sb / (dx[2]*dx[2])) * b[2](i,j,k);
    sten[6] = -(sb / (dx[2]*dx[2])) * b[2](i,j,k+1);
    sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4] + sten[5] + sten[6]);
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

    // zlo
    if (k == boxlo.z) {
        int cdir = Orientation(Direction::z, Orientation::low);
        if (msk[cdir](i,j,k-1) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[2](i,j,k);
            sten[5] = Real(0.0);
            sten[6] += bf2 * b[2](i,j,k);
        }
    }

    // zhi
    if (k == boxhi.z) {
        int cdir = Orientation(Direction::z, Orientation::high);
        if (msk[cdir](i,j,k+1) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[2](i,j,k+1);
            sten[5] += bf2 * b[2](i,j,k+1);
            sten[6] = Real(0.0);
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
        sten[5] = -(sb / (dx[2]*dx[2])) * b[2](i,j,k);
        sten[6] = -(sb / (dx[2]*dx[2])) * b[2](i,j,k+1);
        sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4] + sten[5] + sten[6]);
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

        // zlo
        if (cell_id(i,j,k-1) < 0) {
            int cdir = Orientation(Direction::z, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[2](i,j,k);
            sten[5] = Real(0.0);
            sten[6] += bf2 * b[2](i,j,k);
        }

        // zhi
        if (cell_id(i,j,k+1) < 0) {
            int cdir = Orientation(Direction::z, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[2](i,j,k+1);
            sten[5] += bf2 * b[2](i,j,k+1);
            sten[6] = Real(0.0);
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
void habec_cols (GpuArray<Int,2*AMREX_SPACEDIM+1>& sten, int i, int j, int k,
                 Array4<Int const> const& cell_id)
{
    sten[0] = cell_id(i  ,j  ,k  );
    sten[1] = cell_id(i-1,j  ,k  );
    sten[2] = cell_id(i+1,j  ,k  );
    sten[3] = cell_id(i  ,j-1,k  );
    sten[4] = cell_id(i  ,j+1,k  );
    sten[5] = cell_id(i  ,j  ,k-1);
    sten[6] = cell_id(i  ,j  ,k+1);
}

#ifdef AMREX_USE_EB

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_cols_eb (GpuArray<Int,27>& sten, int i, int j, int k,
                    Array4<Int const> const& cell_id)
{
    sten[0 ] = cell_id(i-1,j-1,k-1);
    sten[1 ] = cell_id(i  ,j-1,k-1);
    sten[2 ] = cell_id(i+1,j-1,k-1);
    sten[3 ] = cell_id(i-1,j  ,k-1);
    sten[4 ] = cell_id(i  ,j  ,k-1);
    sten[5 ] = cell_id(i+1,j  ,k-1);
    sten[6 ] = cell_id(i-1,j+1,k-1);
    sten[7 ] = cell_id(i  ,j+1,k-1);
    sten[8 ] = cell_id(i+1,j+1,k-1);
    sten[9 ] = cell_id(i-1,j-1,k  );
    sten[10] = cell_id(i  ,j-1,k  );
    sten[11] = cell_id(i+1,j-1,k  );
    sten[12] = cell_id(i-1,j  ,k  );
    sten[13] = cell_id(i  ,j  ,k  );
    sten[14] = cell_id(i+1,j  ,k  );
    sten[15] = cell_id(i-1,j+1,k  );
    sten[16] = cell_id(i  ,j+1,k  );
    sten[17] = cell_id(i+1,j+1,k  );
    sten[18] = cell_id(i-1,j-1,k+1);
    sten[19] = cell_id(i  ,j-1,k+1);
    sten[20] = cell_id(i+1,j-1,k+1);
    sten[21] = cell_id(i-1,j  ,k+1);
    sten[22] = cell_id(i  ,j  ,k+1);
    sten[23] = cell_id(i+1,j  ,k+1);
    sten[24] = cell_id(i-1,j+1,k+1);
    sten[25] = cell_id(i  ,j+1,k+1);
    sten[26] = cell_id(i+1,j+1,k+1);
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_ijmat_eb (GpuArray<Real,27>& sten, Array4<Int> const& ncols,
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
                     Array4<Real const> const& apz,
                     Array4<Real const> const& fcx,
                     Array4<Real const> const& fcy,
                     Array4<Real const> const& fcz,
                     Array4<Real const> const& ba,
                     Array4<Real const> const& bcen,
                     Array4<Real const> const& beb)
{
    for (int m = 0; m < 27; ++m) {
        sten[m] = Real(0.0);
    }

    auto& mat_tmp = reinterpret_cast<Array3D<Real,-1,1,-1,1,-1,1>&>(sten);

    GpuArray<Real,AMREX_SPACEDIM> fac{sb/(dx[0]*dx[0]), sb/(dx[1]*dx[1]), sb/(dx[2]*dx[2])};

    if (flag(i,j,k).isRegular())
    {
        mat_tmp(0,0,0) = fac[0]*(b[0](i,j,k)+b[0](i+1,j,k))
                       + fac[1]*(b[1](i,j,k)+b[1](i,j+1,k))
                       + fac[2]*(b[2](i,j,k)+b[2](i,j,k+1));
        mat_tmp(-1, 0, 0) = -fac[0]*b[0](i,j,k);
        mat_tmp( 1, 0, 0) = -fac[0]*b[0](i+1,j,k);
        mat_tmp( 0,-1, 0) = -fac[1]*b[1](i,j,k);
        mat_tmp( 0, 1, 0) = -fac[1]*b[1](i,j+1,k);
        mat_tmp( 0, 0,-1) = -fac[2]*b[2](i,j,k);
        mat_tmp( 0, 0, 1) = -fac[2]*b[2](i,j,k+1);

        if (cell_id(i-1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[0](i,j,k);
            mat_tmp(-1,0,0) = Real(0.0);
            mat_tmp(1,0,0) += bf2*b[0](i,j,k);
        }

        if (cell_id(i+1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[0](i+1,j,k);
            mat_tmp(1,0,0) = Real(0.0);
            mat_tmp(-1,0,0) += bf2*b[0](i+1,j,k);
        }

        if (cell_id(i,j-1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[1](i,j,k);
            mat_tmp(0,-1,0) = Real(0.0);
            mat_tmp(0,1,0) += bf2*b[1](i,j,k);
        }

        if (cell_id(i,j+1,k) < 0) {
            int cdir = Orientation(Direction::y, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[1], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[1](i,j+1,k);
            mat_tmp(0,1,0) = Real(0.0);
            mat_tmp(0,-1,0) += bf2*b[1](i,j+1,k);
        }

        if (cell_id(i,j,k-1) < 0) {
            int cdir = Orientation(Direction::z, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[2](i,j,k);
            mat_tmp(0,0,-1) = Real(0.0);
            mat_tmp(0,0,1) += bf2*b[2](i,j,k);
        }

        if (cell_id(i,j,k+1) < 0) {
            int cdir = Orientation(Direction::z, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[2], bctype[cdir], bcl[cdir], bho);
            mat_tmp(0,0,0) += bf1*b[2](i,j,k+1);
            mat_tmp(0,0,1) = Real(0.0);
            mat_tmp(0,0,-1) += bf2*b[2](i,j,k+1);
        }

        if (sa != Real(0.0)) {
            mat_tmp(0,0,0) += sa*a(i,j,k);
        }
    }
    else if (flag(i,j,k).isSingleValued())
    {
        int cdir = Orientation(Direction::x, Orientation::low);
        Real area = apx(i,j,k);

        if (area > Real(0.0))
        {
            int joff, koff, jj, kk;
            Real fracy, fracz;
            Real f = fac[0];

            if (area != Real(1.0)) {
                joff = static_cast<int>(std::copysign(Real(1.0), fcx(i,j,k,0)));
                koff = static_cast<int>(std::copysign(Real(1.0), fcx(i,j,k,1)));
                jj = j+joff;
                kk = k+koff;
                if (cell_id(i-1,jj,k) < 0 && cell_id(i,jj,k) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcx(i,j,k,0));
                }
                if (cell_id(i-1,j,kk) < 0 && cell_id(i,j,kk) < 0) {
                    fracz = Real(0.0);
                } else {
                    fracz = std::abs(fcx(i,j,k,1));
                }
                if (cell_id(i-1,jj,kk) < 0 && cell_id(i,jj,kk) < 0 && (fracy*fracz) > Real(0.0)) {
                    fracy = Real(0.0);
                    fracz = Real(0.0);
                }
            } else {
                joff = 0;
                koff = 0;
                jj = j;
                kk = k;
                fracy = Real(0.0);
                fracz = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[0], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracy)*(1.0-fracz)*area*b[0](i,j,k);
            // cell(i-1,j,k) is not covered because area > 0
            if (cell_id(i-1,j,k) >= 0) {
                mat_tmp( 0,0,0) += tmp*f;
                mat_tmp(-1,0,0) -= tmp*f;
            } else if (cell_id(i+1,j,k) < 0 || apx(i+1,j,k) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0,0) += tmp*(f+bf1);
                mat_tmp(1,0,0) += tmp*   bf2;
            }

            if (fracy > Real(0.0)) {
                tmp = fracy*(Real(1.0)-fracz)*area*b[0](i,jj,k);
                if (cell_id(i-1,jj,k) >= 0 && cell_id(i,jj,k) >= 0) {
                    mat_tmp(-1,joff,0) -= tmp*f;
                    mat_tmp( 0,joff,0) += tmp*f;
                } else if (cell_id(i+1,jj,k) < 0 || apx(i+1,jj,k) == Real(0.0)) {
                    mat_tmp( 0,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,joff,0) += tmp*(f+bf1);
                    mat_tmp( 1,joff,0) += tmp*   bf2;
                }
            }

            if (fracz > Real(0.0)) {
                tmp = fracz*(Real(1.0)-fracy)*area*b[0](i,j,kk);
                if (cell_id(i-1,j,kk) >= 0 && cell_id(i,j,kk) >= 0) {
                    mat_tmp(-1,0,koff) -= tmp*f;
                    mat_tmp( 0,0,koff) += tmp*f;
                } else if (cell_id(i+1,j,kk) < 0 || apx(i+1,j,kk) == Real(0.0)) {
                    mat_tmp( 0,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,0,koff) += tmp*(f+bf1);
                    mat_tmp( 1,0,koff) += tmp*   bf2;
                }
            }

            if (fracy > Real(0.0) && fracz > Real(0.0)) {
                tmp = fracy*fracz*area*b[0](i,jj,kk);
                if (cell_id(i-1,jj,kk) >= 0 && cell_id(i,jj,kk) >= 0) {
                    mat_tmp(-1,joff,koff) -= tmp*f;
                    mat_tmp( 0,joff,koff) += tmp*f;
                } else if (cell_id(i+1,jj,kk) < 0 || apx(i+1,jj,kk) == Real(0.0)) {
                    mat_tmp( 0,joff,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,joff,koff) += tmp*(f+bf1);
                    mat_tmp( 1,joff,koff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::x, Orientation::high);
        area = apx(i+1,j,k);

        if (area > Real(0.0))
        {
            int joff, koff, jj, kk;
            Real fracy, fracz;
            Real f = fac[0];

            if (area != Real(1.0)) {
                joff = static_cast<int>(std::copysign(Real(1.0), fcx(i+1,j,k,0)));
                koff = static_cast<int>(std::copysign(Real(1.0), fcx(i+1,j,k,1)));
                jj = j+joff;
                kk = k+koff;
                if (cell_id(i,jj,k) < 0 && cell_id(i+1,jj,k) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcx(i+1,j,k,0));
                }
                if (cell_id(i,j,kk) < 0 && cell_id(i+1,j,kk) < 0) {
                    fracz = Real(0.0);
                } else {
                    fracz = std::abs(fcx(i+1,j,k,1));
                }
                if (cell_id(i,jj,kk) < 0 && cell_id(i+1,jj,kk) < 0 && (fracy*fracz) > Real(0.0)) {
                    fracy = Real(0.0);
                    fracz = Real(0.0);
                }
            } else {
                joff = 0;
                koff = 0;
                jj = j;
                kk = k;
                fracy = Real(0.0);
                fracz = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[0], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracy)*(Real(1.0)-fracz)*area*b[0](i+1,j,k);
            if (cell_id(i+1,j,k) >= 0) {
                mat_tmp(0,0,0) += tmp*f;
                mat_tmp(1,0,0) -= tmp*f;
            } else if (cell_id(i-1,j,k) < 0 || apx(i-1,j,k) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp( 0,0,0) += tmp*(f+bf1);
                mat_tmp(-1,0,0) += tmp*   bf2;
            }

            if (fracy > Real(0.0)) {
                tmp = fracy*(Real(1.0)-fracz)*area*b[0](i+1,jj,k);
                if (cell_id(i,jj,k) >= 0 && cell_id(i+1,jj,k) >= 0) {
                    mat_tmp(0,joff,0) += tmp*f;
                    mat_tmp(1,joff,0) -= tmp*f;
                } else if (cell_id(i-1,jj,k) < 0 || apx(i-1,jj,k) == Real(0.0)) {
                    mat_tmp(0,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,joff,0) += tmp*(f+bf1);
                    mat_tmp(-1,joff,0) += tmp*   bf2;
                }
            }

            if (fracz > Real(0.0)) {
                tmp = fracz*(Real(1.0)-fracy)*area*b[0](i+1,j,kk);
                if (cell_id(i,j,kk) >= 0 && cell_id(i+1,j,kk) >= 0) {
                    mat_tmp(0,0,koff) += tmp*f;
                    mat_tmp(1,0,koff) -= tmp*f;
                } else if (cell_id(i-1,j,kk) < 0 || apx(i-1,j,kk) == Real(0.0)) {
                    mat_tmp(0,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,0,koff) += tmp*(f+bf1);
                    mat_tmp(-1,0,koff) += tmp*   bf2;
                }
            }

            if (fracy > Real(0.0) && fracz > Real(0.0)) {
                tmp = fracy*fracz*area*b[0](i+1,jj,kk);
                if (cell_id(i,jj,kk) >= 0 && cell_id(i+1,jj,kk) >= 0) {
                    mat_tmp(0,joff,koff) += tmp*f;
                    mat_tmp(1,joff,koff) -= tmp*f;
                } else if (cell_id(i-1,jj,kk) < 0 || apx(i-1,jj,kk) == Real(0.0)) {
                    mat_tmp(0,joff,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp( 0,joff,koff) += tmp*(f+bf1);
                    mat_tmp(-1,joff,koff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::y, Orientation::low);
        area = apy(i,j,k);

        if (area > Real(0.0))
        {
            int ioff, koff, ii, kk;
            Real fracx, fracz;
            Real f = fac[1];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j,k,0)));
                koff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j,k,1)));
                ii = i+ioff;
                kk = k+koff;
                if (cell_id(ii,j-1,k) < 0 && cell_id(ii,j,k) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcy(i,j,k,0));
                }
                if (cell_id(i,j-1,kk) < 0 && cell_id(i,j,kk) < 0) {
                    fracz = Real(0.0);
                } else {
                    fracz = std::abs(fcy(i,j,k,1));
                }
                if (cell_id(ii,j-1,kk) < 0 && cell_id(ii,j,kk) < 0 && fracx*fracz > Real(0.0)) {
                    fracx = Real(0.0);
                    fracz = Real(0.0);
                }
            } else {
                ioff = 0;
                koff = 0;
                ii = i;
                kk = k;
                fracx = Real(0.0);
                fracz = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[1], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*area*b[1](i,j,k);
            if (cell_id(i,j-1,k) >= 0) {
                mat_tmp(0, 0,0) += tmp*f;
                mat_tmp(0,-1,0) -= tmp*f;
            } else if (cell_id(i,j+1,k) < 0 || apy(i,j+1,k) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0,0) += tmp*(f+bf1);
                mat_tmp(0,1,0) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                tmp = fracx*(Real(1.0)-fracz)*area*b[1](ii,j,k);
                if (cell_id(ii,j-1,k) >= 0 && cell_id(ii,j,k) >= 0) {
                    mat_tmp(ioff,-1,0) -= tmp*f;
                    mat_tmp(ioff, 0,0) += tmp*f;
                } else if (cell_id(ii,j+1,k) < 0 || apy(ii,j+1,k) == Real(0.0)) {
                    mat_tmp(ioff,0,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,0,0) += tmp*(f+bf1);
                    mat_tmp(ioff,1,0) += tmp*   bf2;
                }
            }

            if (fracz > Real(0.0)) {
                tmp = fracz*(Real(1.0)-fracx)*area*b[1](i,j,kk);
                if (cell_id(i,j-1,kk) >= 0 && cell_id(i,j,kk) >= 0) {
                    mat_tmp(0,-1,koff) -= tmp*f;
                    mat_tmp(0, 0,koff) += tmp*f;
                } else if (cell_id(i,j+1,kk) < 0 || apy(i,j+1,kk) == Real(0.0)) {
                    mat_tmp(0,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp(0,0,koff) += tmp*(f+bf1);
                    mat_tmp(0,1,koff) += tmp*   bf2;
                }
            }

            if (fracx > Real(0.0) && fracz > Real(0.0)) {
                tmp = fracx*fracz*area*b[1](ii,j,kk);
                if (cell_id(ii,j-1,kk) >= 0 && cell_id(ii,j,kk) >= 0) {
                    mat_tmp(ioff,-1,koff) -= tmp*f;
                    mat_tmp(ioff, 0,koff) += tmp*f;
                } else if (cell_id(ii,j+1,kk) < 0 || apy(ii,j+1,kk) == Real(0.0)) {
                    mat_tmp(ioff,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,0,koff) += tmp*(f+bf1);
                    mat_tmp(ioff,1,koff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::y, Orientation::high);
        area = apy(i,j+1,k);

        if (area > Real(0.0))
        {
            int ioff, koff, ii, kk;
            Real fracx, fracz;
            Real f = fac[1];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j+1,k,0)));
                koff = static_cast<int>(std::copysign(Real(1.0), fcy(i,j+1,k,1)));
                ii = i+ioff;
                kk = k+koff;
                if (cell_id(ii,j,k) < 0 && cell_id(ii,j+1,k) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcy(i,j+1,k,0));
                }
                if (cell_id(i,j,kk) < 0 && cell_id(i,j+1,kk) < 0) {
                    fracz = Real(0.0);
                } else {
                    fracz = std::abs(fcy(i,j+1,k,1));
                }
                if (cell_id(ii,j,kk) < 0 && cell_id(ii,j+1,kk) < 0 && fracx*fracz > Real(0.0)) {
                    fracx = Real(0.0);
                    fracz = Real(0.0);
                }
            } else {
                ioff = 0;
                koff = 0;
                ii = i;
                kk = k;
                fracx = Real(0.0);
                fracz = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[1], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*area*b[1](i,j+1,k);
            if (cell_id(i,j+1,k) >= 0) {
                mat_tmp(0,0,0) += tmp*f;
                mat_tmp(0,1,0) -= tmp*f;
            } else if (cell_id(i,j-1,k) < 0 || apy(i,j-1,k) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0, 0,0) += tmp*(f+bf1);
                mat_tmp(0,-1,0) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                tmp = fracx*(Real(1.0)-fracz)*area*b[1](ii,j+1,k);
                if (cell_id(ii,j,k) >= 0 && cell_id(ii,j+1,k) >= 0) {
                    mat_tmp(ioff,0,0) += tmp*f;
                    mat_tmp(ioff,1,0) -= tmp*f;
                } else if (cell_id(ii,j-1,k) < 0 || apy(ii,j-1,k) == Real(0.0)) {
                    mat_tmp(ioff,0,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff, 0,0) += tmp*(f+bf1);
                    mat_tmp(ioff,-1,0) += tmp*   bf2;
                }
            }

            if (fracz > Real(0.0)) {
                tmp = fracz*(Real(1.0)-fracx)*area*b[1](i,j+1,kk);
                if (cell_id(i,j,kk) >= 0 && cell_id(i,j+1,kk) >= 0) {
                    mat_tmp(0,0,koff) += tmp*f;
                    mat_tmp(0,1,koff) -= tmp*f;
                } else if (cell_id(i,j-1,kk) < 0 || apy(i,j-1,kk) == Real(0.0)) {
                    mat_tmp(0,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp(0, 0,koff) += tmp*(f+bf1);
                    mat_tmp(0,-1,koff) += tmp*   bf2;
                }
            }

            if (fracx > Real(0.0) && fracz > Real(0.0)) {
                tmp = fracx*fracz*area*b[1](ii,j+1,kk);
                if (cell_id(ii,j,kk) >= 0 && cell_id(ii,j+1,kk) >= 0) {
                    mat_tmp(ioff,1,koff) -= tmp*f;
                    mat_tmp(ioff,0,koff) += tmp*f;
                } else if (cell_id(ii,j-1,kk) < 0 || apy(ii,j-1,kk) == Real(0.0)) {
                    mat_tmp(ioff,0,koff) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff, 0,koff) += tmp*(f+bf1);
                    mat_tmp(ioff,-1,koff) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::z, Orientation::low);
        area = apz(i,j,k);

        if (area > Real(0.0))
        {
            int ioff, joff, ii, jj;
            Real fracx, fracy;
            Real f = fac[2];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcz(i,j,k,0)));
                joff = static_cast<int>(std::copysign(Real(1.0), fcz(i,j,k,1)));
                ii = i+ioff;
                jj = j+joff;
                if (cell_id(ii,j,k-1) < 0 && cell_id(ii,j,k) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcz(i,j,k,0));
                }
                if (cell_id(i,jj,k-1) < 0 && cell_id(i,jj,k) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcz(i,j,k,1));
                }
                if (cell_id(ii,jj,k-1) < 0 && cell_id(ii,jj,k) < 0 && fracx*fracy > Real(0.0)) {
                    fracx = Real(0.0);
                    fracy = Real(0.0);
                }
            } else {
                ioff = 0;
                joff = 0;
                ii = i;
                jj = j;
                fracx = Real(0.0);
                fracy = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[2], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*area*b[2](i,j,k);
            if (cell_id(i,j,k-1) >= 0) {
                mat_tmp(0,0, 0) += tmp*f;
                mat_tmp(0,0,-1) -= tmp*f;
            } else if (cell_id(i,j,k+1) < 0 || apz(i,j,k+1) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0,0) += tmp*(f+bf1);
                mat_tmp(0,0,1) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                tmp = fracx*(Real(1.0)-fracy)*area*b[2](ii,j,k);
                if (cell_id(ii,j,k-1) >= 0 && cell_id(ii,j,k) >= 0) {
                    mat_tmp(ioff,0,-1) -= tmp*f;
                    mat_tmp(ioff,0, 0) += tmp*f;
                } else if (cell_id(ii,j,k+1) < 0 || apz(ii,j,k+1) == Real(0.0)) {
                    mat_tmp(ioff,0,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,0,0) += tmp*(f+bf1);
                    mat_tmp(ioff,0,1) += tmp*   bf2;
                }
            }

            if (fracy > Real(0.0)) {
                tmp = fracy*(Real(1.0)-fracx)*area*b[2](i,jj,k);
                if (cell_id(i,jj,k-1) >= 0 && cell_id(i,jj,k) >= 0) {
                    mat_tmp(0,joff,-1) -= tmp*f;
                    mat_tmp(0,joff, 0) += tmp*f;
                } else if (cell_id(i,jj,k+1) < 0 || apz(i,jj,k+1) == Real(0.0)) {
                    mat_tmp(0,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(0,joff,0) += tmp*(f+bf1);
                    mat_tmp(0,joff,1) += tmp*   bf2;
                }
            }

            if (fracx > Real(0.0) && fracy > Real(0.0)) {
                tmp = fracx*fracy*area*b[2](ii,jj,k);
                if (cell_id(ii,jj,k-1) >= 0 && cell_id(ii,jj,k) >= 0) {
                    mat_tmp(ioff,joff,-1) -= tmp*f;
                    mat_tmp(ioff,joff, 0) += tmp*f;
                } else if (cell_id(ii,jj,k+1) < 0 || apz(ii,jj,k+1) == Real(0.0)) {
                    mat_tmp(ioff,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,joff,0) += tmp*(f+bf1);
                    mat_tmp(ioff,joff,1) += tmp*   bf2;
                }
            }
        }

        cdir = Orientation(Direction::z, Orientation::high);
        area = apz(i,j,k+1);

        if (area > Real(0.0))
        {
            int ioff, joff, ii, jj;
            Real fracx, fracy;
            Real f = fac[2];

            if (area != Real(1.0)) {
                ioff = static_cast<int>(std::copysign(Real(1.0), fcz(i,j,k+1,0)));
                joff = static_cast<int>(std::copysign(Real(1.0), fcz(i,j,k+1,1)));
                ii = i+ioff;
                jj = j+joff;
                if (cell_id(ii,j,k) < 0 && cell_id(ii,j,k+1) < 0) {
                    fracx = Real(0.0);
                } else {
                    fracx = std::abs(fcz(i,j,k+1,0));
                }
                if (cell_id(i,jj,k) < 0 && cell_id(i,jj,k+1) < 0) {
                    fracy = Real(0.0);
                } else {
                    fracy = std::abs(fcz(i,j,k+1,1));
                }
                if (cell_id(ii,jj,k) < 0 && cell_id(ii,jj,k+1) < 0 && fracx*fracy > Real(0.0)) {
                    fracx = Real(0.0);
                    fracy = Real(0.0);
                }
            } else {
                ioff = 0;
                joff = 0;
                ii = i;
                jj = j;
                fracx = Real(0.0);
                fracy = Real(0.0);
            }

            Real bf1, bf2, bflo;
            detail::comp_bflo(bf1, bf2, bflo, sb, dx[2], bctype[cdir], bcl[cdir], bho);

            Real tmp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*area*b[2](i,j,k+1);
            if (cell_id(i,j,k+1) >= 0) {
                mat_tmp(0,0,0) += tmp*f;
                mat_tmp(0,0,1) -= tmp*f;
            } else if (cell_id(i,j,k-1) < 0 || apz(i,j,k-1) == Real(0.0)) {
                mat_tmp(0,0,0) += tmp*(f+bflo);
            } else {
                mat_tmp(0,0, 0) += tmp*(f+bf1);
                mat_tmp(0,0,-1) += tmp*   bf2;
            }

            if (fracx > Real(0.0)) {
                tmp = fracx*(Real(1.0)-fracy)*area*b[2](ii,j,k+1);
                if (cell_id(ii,j,k) >= 0 && cell_id(ii,j,k+1) >= 0) {
                    mat_tmp(ioff,0,0) += tmp*f;
                    mat_tmp(ioff,0,1) -= tmp*f;
                } else if (cell_id(ii,j,k-1) < 0 || apz(ii,j,k-1) == Real(0.0)) {
                    mat_tmp(ioff,0,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,0, 0) += tmp*(f+bf1);
                    mat_tmp(ioff,0,-1) += tmp*   bf2;
                }
            }

            if (fracy > Real(0.0)) {
                tmp = fracy*(Real(1.0)-fracx)*area*b[2](i,jj,k+1);
                if (cell_id(i,jj,k) >= 0 && cell_id(i,jj,k+1) >= 0) {
                    mat_tmp(0,joff,0) += tmp*f;
                    mat_tmp(0,joff,1) -= tmp*f;
                } else if (cell_id(i,jj,k-1) < 0 || apz(i,jj,k-1) == Real(0.0)) {
                    mat_tmp(0,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(0,joff, 0) += tmp*(f+bf1);
                    mat_tmp(0,joff,-1) += tmp*   bf2;
                }
            }

            if (fracx > Real(0.0) && fracy > Real(0.0)) {
                tmp = fracx*fracy*area*b[2](ii,jj,k+1);
                if (cell_id(ii,jj,k+1) >= 0 && cell_id(ii,jj,k) >= 0) {
                    mat_tmp(ioff,joff,1) -= tmp*f;
                    mat_tmp(ioff,joff,0) += tmp*f;
                } else if (cell_id(ii,jj,k-1) < 0 || apz(ii,jj,k-1) == Real(0.0)) {
                    mat_tmp(ioff,joff,0) += tmp*(f+bflo);
                } else {
                    mat_tmp(ioff,joff, 0) += tmp*(f+bf1);
                    mat_tmp(ioff,joff,-1) += tmp*   bf2;
                }
            }
        }

        if (beb) {
            Real anorm = std::sqrt((apx(i,j,k) - apx(i+1,j,k))*(apx(i,j,k) - apx(i+1,j,k))
                                 + (apy(i,j,k) - apy(i,j+1,k))*(apy(i,j,k) - apy(i,j+1,k))
                                 + (apz(i,j,k) - apz(i,j,k+1))*(apz(i,j,k) - apz(i,j,k+1)));

            Real anorminv = Real(1.0)/anorm;
            Real anrmx = (apx(i,j,k) - apx(i+1,j,k))*anorminv;
            Real anrmy = (apy(i,j,k) - apy(i,j+1,k))*anorminv;
            Real anrmz = (apz(i,j,k) - apz(i,j,k+1))*anorminv;
            Real sx   = std::copysign(Real(1.0),anrmx);
            Real sy   = std::copysign(Real(1.0),anrmy);
            Real sz   = std::copysign(Real(1.0),anrmz);
            Real bctx = bcen(i,j,k,0);
            Real bcty = bcen(i,j,k,1);
            Real bctz = bcen(i,j,k,2);
            Real dg   = get_dx_eb(vfrc(i,j,k)) / amrex::max(std::abs(anrmx), std::abs(anrmy), std::abs(anrmz));
            Real gx   = bctx - dg*anrmx;
            Real gy   = bcty - dg*anrmy;
            Real gz   = bctz - dg*anrmz;
            int ioff = -static_cast<int>(sx);
            int joff = -static_cast<int>(sy);
            int koff = -static_cast<int>(sz);
            gx *= sx;
            gy *= sy;
            gz *= sz;
            Real gxy  = gx*gy;
            Real gxz  = gx*gz;
            Real gyz  = gy*gz;
            Real gxyz = gx*gy*gz;
            Array1D<Real,0,7> phig1{Real(1.0) + gx + gy + gz + gxy + gxz + gyz + gxyz,
                                              - gx           - gxy - gxz       - gxyz,
                                                   - gy      - gxy       - gyz - gxyz,
                                                        - gz       - gxz - gyz - gxyz,
                                                             + gxy             + gxyz,
                                                                   + gxz       + gxyz,
                                                                         + gyz + gxyz,
                                                                               - gxyz};
            Array1D<Real,0,7> feb;
            for (int ii=0; ii<8; ii++) {
                feb(ii) = -phig1(ii) * (ba(i,j,k) * beb(i,j,k) / dg);
            }
            mat_tmp(0   , 0  , 0  ) -= feb(0)*fac[0];
            mat_tmp(ioff, 0  , 0  ) -= feb(1)*fac[0];
            mat_tmp(0   ,joff, 0  ) -= feb(2)*fac[0];
            mat_tmp(0   , 0  ,koff) -= feb(3)*fac[0];
            mat_tmp(ioff,joff, 0  ) -= feb(4)*fac[0];
            mat_tmp(ioff, 0  ,koff) -= feb(5)*fac[0];
            mat_tmp(0   ,joff,koff) -= feb(6)*fac[0];
            mat_tmp(ioff,joff,koff) -= feb(7)*fac[0];
        }

        for (int kk=-1; kk<=1; kk++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            mat_tmp(ii,jj,kk) *= Real(1.0)/vfrc(i,j,k);
        }}}

        if (sa != Real(0.0)) {
            mat_tmp(0,0,0) += sa*a(i,j,k);
        }
    }
    else
    {
        for (int kk=-1; kk<=1; kk++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            mat_tmp(ii,jj,kk) = Real(0.0);
        }}}
        mat_tmp(0,0,0) = Real(1.0);
    }

    diaginv(i,j,k) = Real(1.0)/mat_tmp(0,0,0);
    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        mat_tmp(ii,jj,kk) *= diaginv(i,j,k);
    }}}
    mat_tmp(0,0,0) = Real(1.0);

    ncols(i,j,k) = 0;
    for (int m = 0; m < 27; ++m) {
        ncols(i,j,k) += (sten[m] != Real(0.0));
    }
}

#endif

}
#endif
