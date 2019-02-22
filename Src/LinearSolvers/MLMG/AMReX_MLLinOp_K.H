#ifndef AMREX_MLLINOP_K_H_
#define AMREX_MLLINOP_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BoundCond.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LOUtil_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mllinop_apply_bc_x (Box const& box, int blen, Array4<Real> const& phi,
                         Array4<int const> const& masklo, Array4<int const> const& maskhi,
                         BoundCond bctlo, BoundCond bcthi, Real bcllo, Real bclhi,
                         Array4<Real const> const& bcvallo, Array4<Real const> const& bcvalhi,
                         int maxorder, Real dxinv, int inhomog, int ncomp)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    int i = lo.x;
    switch (bctlo) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i+1,j,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i+1,j,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.x-lo.x+2, maxorder);
        x[0] = -bcllo * dxinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (masklo(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i+m,j,k,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvallo(i,j,k,n)*coef[0];
                        }
                    }
                }
            }
        }
        break;
    }
    default: {}
    }

    i = lo.x + blen + 1;
    switch (bcthi) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i-1,j,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i-1,j,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.x-lo.x+2, maxorder);
        x[0] = -bclhi * dxinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    if (maskhi(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i-m,j,k,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvalhi(i,j,k,n)*coef[0];
                        }
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
void mllinop_apply_bc_y (Box const& box, int blen, Array4<Real> const& phi,
                         Array4<int const> const& masklo, Array4<int const> const& maskhi,
                         BoundCond bctlo, BoundCond bcthi, Real bcllo, Real bclhi,
                         Array4<Real const> const& bcvallo, Array4<Real const> const& bcvalhi,
                         int maxorder, Real dyinv, int inhomog, int ncomp)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    int j = lo.y;
    switch (bctlo) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i,j+1,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i,j+1,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.y-lo.y+2, maxorder);
        x[0] = -bcllo * dyinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i,j+m,k,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvallo(i,j,k,n)*coef[0];
                        }
                    }
                }
            }
        }
        break;
    }
    default: {}
    }

    j = lo.y + blen + 1;
    switch (bcthi) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i,j-1,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i,j-1,k,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.y-lo.y+2, maxorder);
        x[0] = -bclhi * dyinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i,j-m,k,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvalhi(i,j,k,n)*coef[0];
                        }
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
void mllinop_apply_bc_z (Box const& box, int blen, Array4<Real> const& phi,
                         Array4<int const> const& masklo, Array4<int const> const& maskhi,
                         BoundCond bctlo, BoundCond bcthi, Real bcllo, Real bclhi,
                         Array4<Real const> const& bcvallo, Array4<Real const> const& bcvalhi,
                         int maxorder, Real dzinv, int inhomog, int ncomp)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    int k = lo.z;
    switch (bctlo) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i,j,k+1,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i,j,k+1,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.z-lo.z+2, maxorder);
        x[0] = -bcllo * dzinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (masklo(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i,j,k+m,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvallo(i,j,k,n)*coef[0];
                        }
                    }
                }
            }
        }
        break;
    }
    default: {}
    }

    k = lo.z + blen + 1;
    switch (bcthi) {
    case AMREX_LO_NEUMANN:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = phi(i,j,k-1,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_REFLECT_ODD:
    {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        phi(i,j,k,n) = -phi(i,j,k-1,n);
                    }
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        Real x[4], coef[4];
        const int NX = amrex::min(hi.x-lo.x+2, maxorder);
        x[0] = -bclhi * dzinv;
        for (int m = 1; m < NX; ++m) {
            x[m] = m - 0.5;
        }
        poly_interp_coeff(-0.5, x, NX, coef);
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (maskhi(i,j,k) > 0) {
                        Real s = 0.0;
                        for (int m = 1; m < NX; ++m) {
                            s += phi(i,j,k-m,n) * coef[m];
                        }
                        phi(i,j,k,n) = s;
                        if (inhomog) {
                            phi(i,j,k,n) += bcvalhi(i,j,k,n)*coef[0];
                        }
                    }
                }
            }
        }
        break;
    }
    default: {}
    }
}

}

#endif
