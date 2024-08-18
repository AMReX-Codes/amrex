#ifndef AMREX_HYPRE_ML_ABECLAP_K_H_
#define AMREX_HYPRE_ML_ABECLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LOUtil_K.H>
#include <AMReX_REAL.H>

#include <HYPRE_utilities.h>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_mat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, int i, int j, int k,
                       Dim3 const& boxlo, Dim3 const& boxhi,
                       Real sa, Array4<Real const> const& a,
                       Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                       GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                       GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                       GpuArray<Real,AMREX_SPACEDIM*2> const& bcl,
                       GpuArray<Array4<int const>, AMREX_SPACEDIM*2> const& bcmsk,
                       GpuArray<Array4<Real const>, AMREX_SPACEDIM*2> const& bcval,
                       GpuArray<Array4<Real>, AMREX_SPACEDIM*2> const& bcrhs,
                       int level, IntVect const& fixed_pt)
{
    Real bxm = b[0] ? b[0](i  ,j  ,k  ) : Real(1.0);
    Real bxp = b[0] ? b[0](i+1,j  ,k  ) : Real(1.0);
    Real bym = b[1] ? b[1](i  ,j  ,k  ) : Real(1.0);
    Real byp = b[1] ? b[1](i  ,j+1,k  ) : Real(1.0);
#if (AMREX_SPACEDIM > 2)
    Real bzm = b[2] ? b[2](i  ,j  ,k  ) : Real(1.0);
    Real bzp = b[2] ? b[2](i  ,j  ,k+1) : Real(1.0);
#endif
    Real ac = a ? a(i,j,k) : Real(0.0);

    sten[1] = -(sb / (dx[0]*dx[0])) * bxm;
    sten[2] = -(sb / (dx[0]*dx[0])) * bxp;
    sten[3] = -(sb / (dx[1]*dx[1])) * bym;
    sten[4] = -(sb / (dx[1]*dx[1])) * byp;
#if (AMREX_SPACEDIM == 2)
    sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4]) + sa*ac;
#else
    sten[5] = -(sb / (dx[2]*dx[2])) * bzm;
    sten[6] = -(sb / (dx[2]*dx[2])) * bzp;
    sten[0] = -(sten[1] + sten[2] + sten[3] + sten[4] + sten[5] + sten[6]) + sa*ac;
#endif

    // xlo
    if (i == boxlo.x) {
        int const cdir = Orientation(Direction::x, Orientation::low);
        int const bcmk = bcmsk[cdir](i-1,j,k);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[0]*Real(0.5), dx[0]*Real(1.5)};
                poly_interp_coeff<3>(dx[0]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[0]*dx[0])) * bxm;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i-1,j,k);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[1] = Real(0.0);
            sten[2] -= fac * cc[2];
        }
    }

    // xhi
    if (i == boxhi.x) {
        int const cdir = Orientation(Direction::x, Orientation::high);
        int const bcmk = bcmsk[cdir](i+1,j,k);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[0]*Real(0.5), dx[0]*Real(1.5)};
                poly_interp_coeff<3>(dx[0]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[0]*dx[0])) * bxp;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i+1,j,k);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[1] -= fac * cc[2];
            sten[2] = Real(0.0);
        }
    }

    // ylo
    if (j == boxlo.y) {
        int const cdir = Orientation(Direction::y, Orientation::low);
        int const bcmk = bcmsk[cdir](i,j-1,k);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[1]*Real(0.5), dx[1]*Real(1.5)};
                poly_interp_coeff<3>(dx[1]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[1]*dx[1])) * bym;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i,j-1,k);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[3] = Real(0.0);
            sten[4] -= fac * cc[2];
        }
    }

    // yhi
    if (j == boxhi.y) {
        int const cdir = Orientation(Direction::y, Orientation::high);
        int const bcmk = bcmsk[cdir](i,j+1,k);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[1]*Real(0.5), dx[1]*Real(1.5)};
                poly_interp_coeff<3>(dx[1]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[1]*dx[1])) * byp;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i,j+1,k);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[3] -= fac * cc[2];
            sten[4] = Real(0.0);
        }
    }

#if (AMREX_SPACEDIM > 2)

    // zlo
    if (k == boxlo.z) {
        int const cdir = Orientation(Direction::z, Orientation::low);
        int const bcmk = bcmsk[cdir](i,j,k-1);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[2]*Real(0.5), dx[2]*Real(1.5)};
                poly_interp_coeff<3>(dx[2]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[2]*dx[2])) * bzm;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i,j,k-1);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[5] = Real(0.0);
            sten[6] -= fac * cc[2];
        }
    }

    // zhi
    if (k == boxhi.z) {
        int const cdir = Orientation(Direction::z, Orientation::high);
        int const bcmk = bcmsk[cdir](i,j,k+1);
        if (bcmk > 0) {
            int const bct = bctype[cdir];
            Real cc[3];
            if (bct == AMREX_LO_DIRICHLET) {
                Real xx[3] = {-bcl[cdir], dx[2]*Real(0.5), dx[2]*Real(1.5)};
                poly_interp_coeff<3>(dx[2]*Real(-0.5), xx, cc);
            } else { // Neumann
                cc[0] = Real(0.0);
                cc[1] = Real(1.0);
                cc[2] = Real(0.0);
            }
            Real fac = (sb / (dx[2]*dx[2])) * bzp;
            if (bct == AMREX_LO_DIRICHLET && (level == 0 || bcmk == 2)) {
                // bcmk == 2 means outside the domain.
                // We need to modify RHS at external Dirichlet boundaries.
                bcrhs[cdir](i,j,k) = fac * cc[0] * bcval[cdir](i,j,k+1);
            } else {
                bcrhs[cdir](i,j,k) = Real(0.0);
            }
            sten[0] -= fac * cc[1];
            sten[5] -= fac * cc[2];
            sten[6] = Real(0.0);
        }
    }

#endif

    if (fixed_pt == IntVect(AMREX_D_DECL(i,j,k))) {
        for (int n = 1; n < 2*AMREX_SPACEDIM+1; ++n) {
            sten[n] = Real(0.0);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_rhs (int i, int j, int k, Dim3 const& boxlo, Dim3 const& boxhi,
                       Array4<Real> const& rhs1,
                       Array4<Real const> const& rhs0,
                       GpuArray<Array4<int const>, AMREX_SPACEDIM*2> const& bcmsk,
                       GpuArray<Array4<Real const>, AMREX_SPACEDIM*2> const& bcrhs)
{
    rhs1(i,j,k) = rhs0(i,j,k);

    // xlo
    if (i == boxlo.x) {
        int cdir = Orientation(Direction::x, Orientation::low);
        if (bcmsk[cdir](i-1,j,k) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

    // xhi
    if (i == boxhi.x) {
        int cdir = Orientation(Direction::x, Orientation::high);
        if (bcmsk[cdir](i+1,j,k) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

    // ylo
    if (j == boxlo.y) {
        int cdir = Orientation(Direction::y, Orientation::low);
        if (bcmsk[cdir](i,j-1,k) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

    // yhi
    if (j == boxhi.y) {
        int cdir = Orientation(Direction::y, Orientation::high);
        if (bcmsk[cdir](i,j+1,k) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

#if (AMREX_SPACEDIM > 2)

    // zlo
    if (k == boxlo.z) {
        int cdir = Orientation(Direction::z, Orientation::low);
        if (bcmsk[cdir](i,j,k-1) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

    // zhi
    if (k == boxhi.z) {
        int cdir = Orientation(Direction::z, Orientation::high);
        if (bcmsk[cdir](i,j,k+1) > 0) {
            rhs1(i,j,k) += bcrhs[cdir](i,j,k);
        }
    }

#endif
}

}

#if (AMREX_SPACEDIM == 2)
#include <AMReX_HypreMLABecLap_2D_K.H>
#else
#include <AMReX_HypreMLABecLap_3D_K.H>
#endif

#endif
