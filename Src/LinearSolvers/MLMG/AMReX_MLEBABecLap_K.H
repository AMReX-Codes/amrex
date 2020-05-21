#ifndef AMREX_MLEBABECLAP_K_H_
#define AMREX_MLEBABECLAP_K_H_

#include <AMReX_MLLinOp_K.H>

#include <AMReX_EBCellFlag.H>

namespace amrex { namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real get_dx_eb (Real kappa) noexcept {
        return amrex::max(0.3_rt,(kappa*kappa-0.25_rt)/(2.0_rt*kappa));
    }
}}

#if (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBABecLap_2D_K.H>
#else
#include <AMReX_MLEBABecLap_3D_K.H>
#endif

namespace amrex {

// note that the mask in these functions is different from masks in bndry registers
// 1 means valid data, 0 means invalid data

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_apply_bc_x (int side, Box const& box, int blen,
                             Array4<Real> const& phi,
                             Array4<int const> const& mask,
                             Array4<Real const> const& area,
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
                if (mask(i,j,k) == 0 and mask(i+s,j,k) == 1) {
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
                if (mask(i,j,k) == 0 and mask(i+s,j,k) == 1) {
                    phi(i,j,k,icomp) = -phi(i+s,j,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        const int NX = amrex::min(blen+1, maxorder);
        GpuArray<Real,4> x{-bcl * dxinv, 0.5_rt, 1.5_rt, 2.5_rt};
        Array2D<Real, 0, 3, 0, 2> coef{};
        for (int r = 0; r <= maxorder-2; ++r) {
            poly_interp_coeff(-0.5_rt, &x[0], NX, &(coef(0,r)));
        }
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                if (mask(i,j,k) == 0 and mask(i+s,j,k) == 1) {
                    int order = 1;
                    for (int r = 0; r <= NX-2; ++r) {
                        if (area(i+(1-side)+s*r,j,k) > 0.0_rt) {
                            ++order;
                        } else {
                            break;
                        }
                    }
                    if (order == 1) {
                        if (inhomog) {
                            phi(i,j,k,icomp) = bcval(i,j,k,icomp);
                        } else {
                            phi(i,j,k,icomp) = 0.0_rt;
                        }
                    } else {
                        Real tmp = 0.0_rt;
                        for (int m = 1; m < order; ++m) {
                            tmp += phi(i+m*s,j,k,icomp) * coef(m,order-2);
                        }
                        phi(i,j,k,icomp) = tmp;
                        if (inhomog) {
                            phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef(0,order-2);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_apply_bc_y (int side, Box const& box, int blen,
                             Array4<Real> const& phi,
                             Array4<int const> const& mask,
                             Array4<Real const> const& area,
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
                if (mask(i,j,k) == 0 and mask(i,j+s,k) == 1) {
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
                if (mask(i,j,k) == 0 and mask(i,j+s,k) == 1) {
                    phi(i,j,k,icomp) = -phi(i,j+s,k,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        const int NX = amrex::min(blen+1, maxorder);
        GpuArray<Real,4> x{-bcl * dyinv, 0.5_rt, 1.5_rt, 2.5_rt};
        Array2D<Real, 0, 3, 0, 2> coef{};
        for (int r = 0; r <= maxorder-2; ++r) {
            poly_interp_coeff(-0.5_rt, &x[0], NX, &(coef(0,r)));
        }
        for     (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) == 0 and mask(i,j+s,k) == 1) {
                    int order = 1;
                    for (int r = 0; r <= NX-2; ++r) {
                        if (area(i,j+(1-side)+s*r,k) > 0.0_rt) {
                            ++order;
                        } else {
                            break;
                        }
                    }
                    if (order == 1) {
                        if (inhomog) {
                            phi(i,j,k,icomp) = bcval(i,j,k,icomp);
                        } else {
                            phi(i,j,k,icomp) = 0.0_rt;
                        }
                    } else {
                        Real tmp = 0.0_rt;
                        for (int m = 1; m < order; ++m) {
                            tmp += phi(i,j+m*s,k,icomp) * coef(m,order-2);
                        }
                        phi(i,j,k,icomp) = tmp;
                        if (inhomog) {
                            phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef(0,order-2);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_apply_bc_z (int side, Box const& box, int blen,
                             Array4<Real> const& phi,
                             Array4<int const> const& mask,
                             Array4<Real const> const& area,
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
                if (mask(i,j,k) == 0 and mask(i,j,k+s) == 1) {
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
                if (mask(i,j,k) == 0 and mask(i,j,k+s) == 1) {
                    phi(i,j,k,icomp) = -phi(i,j,k+s,icomp);
                }
            }
        }
        break;
    }
    case AMREX_LO_DIRICHLET:
    {
        const int NX = amrex::min(blen+1, maxorder);
        GpuArray<Real,4> x{-bcl * dzinv, 0.5_rt, 1.5_rt, 2.5_rt};
        Array2D<Real, 0, 3, 0, 2> coef{};
        for (int r = 0; r <= maxorder-2; ++r) {
            poly_interp_coeff(-0.5_rt, &x[0], NX, &(coef(0,r)));
        }
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mask(i,j,k) == 0 and mask(i,j,k+s) == 1) {
                    int order = 1;
                    for (int r = 0; r <= NX-2; ++r) {
                        if (area(i,j,k+(1-side)+s*r) > 0.0_rt) {
                            ++order;
                        } else {
                            break;
                        }
                    }
                    if (order == 1) {
                        if (inhomog) {
                            phi(i,j,k,icomp) = bcval(i,j,k,icomp);
                        } else {
                            phi(i,j,k,icomp) = 0.0_rt;
                        }
                    } else {
                        Real tmp = 0.0_rt;
                        for (int m = 1; m < order; ++m) {
                            tmp += phi(i,j,k+m*s,icomp) * coef(m,order-2);
                        }
                        phi(i,j,k,icomp) = tmp;
                        if (inhomog) {
                            phi(i,j,k,icomp) += bcval(i,j,k,icomp)*coef(0,order-2);
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

