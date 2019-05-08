#ifndef AMREX_MLLINOP_K_H_
#define AMREX_MLLINOP_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BoundCond.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LOUtil_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_apply_bc_x (int side, Box const& box, int blen,
                         Array4<Real> const& phi,
                         Array4<int const> const& mask,
                         BoundCond bct, Real bcl,
                         Array4<Real const> const& bcval,
                         int maxorder, Real dxinv, int inhomog, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int i = lo.x; // boundary cell
    const int s = 1-2*side;  // +1 for lo and -1 for hi
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = phi(i+s,j,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = -phi(i+s,j,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dxinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                if (mask(i,j,k) > 0) {
                    Real tmp = 0.0;
                    for (int m = 1; m < NX; ++m) {
                        tmp += phi(i+m*s,j,k,icomp) * coef[m];
                    }
                    phi(i,j,k,icomp) = tmp;
                    if (inhomog) {
                        phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef[0];
                    }
                }
            }
        }
        break;
    }
    default: {}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_apply_bc_y (int side, Box const& box, int blen,
                         Array4<Real> const& phi,
                         Array4<int const> const& mask,
                         BoundCond bct, Real bcl,
                         Array4<Real const> const& bcval,
                         int maxorder, Real dyinv, int inhomog, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int j = lo.y; // boundary cell
    const int s = 1-2*side; // +1 for lo and -1 for hi
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = phi(i,j+s,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = -phi(i,j+s,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dyinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    Real tmp = 0.0;
                    for (int m = 1; m < NX; ++m) {
                        tmp += phi(i,j+m*s,k,icomp) * coef[m];
                    }
                    phi(i,j,k,icomp) = tmp;
                    if (inhomog) {
                        phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef[0];
                    }
                }
            }
        }
        break;
    }
    default: {}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_apply_bc_z (int side, Box const& box, int blen,
                         Array4<Real> const& phi,
                         Array4<int const> const& mask,
                         BoundCond bct, Real bcl,
                         Array4<Real const> const& bcval,
                         int maxorder, Real dzinv, int inhomog, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int k = lo.z; // boundary cell
    const int s = 1-2*side; // +1 for lo and -1 for hi
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = phi(i,j,k+s,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    phi(i,j,k,icomp) = -phi(i,j,k+s,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dzinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) > 0) {
                    Real tmp = 0.0;
                    for (int m = 1; m < NX; ++m) {
                        tmp += phi(i,j,k+m*s,icomp) * coef[m];
                    }
                    phi(i,j,k,icomp) = tmp;
                    if (inhomog) {
                        phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef[0];
                    }
                }
            }
        }
        break;
    }
    default: {}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_comp_interp_coef0_x (int side, Box const& box, int blen,
                                  Array4<Real> const& f,
                                  Array4<int const> const& mask,
                                  BoundCond bct, Real bcl,
                                  int maxorder, Real dxinv, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int ib = lo.x; // boundary cell
    const int ii = lo.x + (1-2*side); // interior cell
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                f(ii,j,k,icomp) = 1.0;
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                f(ii,j,k,icomp) = (mask(ib,j,k) > 0) ? 1.0 : 0.0;
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dxinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                f(ii,j,k,icomp) = (mask(ib,j,k) > 0) ? coef[1] : 0.0;
            }
        }
        break;
    }
    default: {}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_comp_interp_coef0_y (int side, Box const& box, int blen,
                                  Array4<Real> const& f,
                                  Array4<int const> const& mask,
                                  BoundCond bct, Real bcl,
                                  int maxorder, Real dyinv, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int jb = lo.y; // boundary cell
    const int ji = lo.y + (1-2*side); // interior cell
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,ji,k,icomp) = 1.0;
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,ji,k,icomp) = (mask(i,jb,k) > 0) ? 1.0 : 0.0;
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dyinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,ji,k,icomp) = (mask(i,jb,k) > 0) ? coef[1] : 0.0;
            }
        }
        break;
    }
    default: {}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_comp_interp_coef0_z (int side, Box const& box, int blen,
                                  Array4<Real> const& f,
                                  Array4<int const> const& mask,
                                  BoundCond bct, Real bcl,
                                  int maxorder, Real dzinv, int icomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const int kb = lo.z; // bound cell
    const int ki = lo.z + (1-2*side); // interior cell
    switch (bct) {
    case AMREX_LO_NEUMANN:
    {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,j,ki,icomp) = 1.0;
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,j,ki,icomp) = (mask(i,j,kb) > 0) ? 1.0 : 0.0;
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(blen+1, maxorder);
        x[0] = -bcl * dzinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                f(i,j,ki,icomp) = (mask(i,j,kb) > 0) ? coef[1] : 0.0;
            }
        }
        break;
    }
    default: {}
    }
}

}

#endif
