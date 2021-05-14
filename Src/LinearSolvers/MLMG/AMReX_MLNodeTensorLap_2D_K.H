#ifndef AMREX_MLNODETENSORLAP_2D_K_H_
#define AMREX_MLNODETENSORLAP_2D_K_H_
#include <AMReX_Config.H>

namespace amrex {

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_x (Array4<Real const> const& crse, int ic, int jc) noexcept
    {
        return (crse(ic,jc,0)+crse(ic+1,jc,0))*Real(0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_y (Array4<Real const> const& crse, int ic, int jc) noexcept
    {
        return (crse(ic,jc,0)+crse(ic,jc+1,0))*Real(0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_face_xy (Array4<Real const> const& crse, int ic, int jc) noexcept
    {
        return (ts_interp_line_y(crse,ic  ,jc  ) +
                ts_interp_line_y(crse,ic+1,jc  ) +
                ts_interp_line_x(crse,ic  ,jc  ) +
                ts_interp_line_x(crse,ic  ,jc+1)) * Real(0.25);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int i, int j, int, Array4<Real> const& fine,
                          Array4<Real const> const& crse, Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += ts_interp_face_xy(crse,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += ts_interp_line_x(crse,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += ts_interp_line_y(crse,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const& b, Array4<Real> const& y, Array4<Real const> const& x,
                      GpuArray<Real,3> const& s,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h11 = dxinv[1]*dxinv[1];
    amrex::LoopConcurrent(b, [=] (int i, int j, int k) noexcept
    {
        y(i,j,k) = x(i-1,j-1,0) * (Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
            +      x(i-1,j  ,0) * (Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2])
            +      x(i-1,j+1,0) * (Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
            +      x(i  ,j-1,0) * (Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2])
            +      x(i  ,j  ,0) * (Real(-4./3.)*h00*s[0] + Real(-4./3.)*h11*s[2])
            +      x(i  ,j+1,0) * (Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2])
            +      x(i+1,j-1,0) * (Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
            +      x(i+1,j  ,0) * (Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2])
            +      x(i+1,j+1,0) * (Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2]);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_gauss_seidel (Box const& b, Array4<Real> const& sol,
                             Array4<Real const> const& rhs, Array4<int const> const& msk,
                             GpuArray<Real,3> const& s,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h11 = dxinv[1]*dxinv[1];
    amrex::Loop(b, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (Real(-4./3.)*h00*s[0] + Real(-4./3.)*h11*s[2]);
            Real Ax = sol(i-1,j-1,0) * (Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
                +     sol(i-1,j  ,0) * (Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2])
                +     sol(i-1,j+1,0) * (Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
                +     sol(i  ,j-1,0) * (Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2])
                +     sol(i  ,j  ,0) * s0
                +     sol(i  ,j+1,0) * (Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2])
                +     sol(i+1,j-1,0) * (Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2])
                +     sol(i+1,j  ,0) * (Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2])
                +     sol(i+1,j+1,0) * (Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2]);
            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_normalize (Box const& b, Array4<Real> const& phi,
                          Array4<int const> const& msk, GpuArray<Real,3> const& s,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h11 = dxinv[1]*dxinv[1];
    amrex::Loop(b, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            Real s0 = (Real(-4./3.)*h00*s[0] + Real(-4./3.)*h11*s[2]);
            phi(i,j,k) /= s0;
        }
    });
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)

template <typename HypreInt, typename AtomicInt>
void mlndtslap_fill_ijmatrix_cpu (Box const& ndbx,
                                  Array4<AtomicInt const> const& gid,
                                  Array4<int const> const& lid,
                                  HypreInt* const ncols, HypreInt* const cols, Real* const mat,
                                  GpuArray<Real,3> const& s,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h11 = dxinv[1]*dxinv[1];

    constexpr auto gidmax = std::numeric_limits<AtomicInt>::max();
    HypreInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            HypreInt nelems_old = nelems;

            cols[nelems] = gid(i,j,k);
            mat[nelems] = Real(-4./3.)*h00*s[0] + Real(-4./3.)*h11*s[2];
            ++nelems;

            if                (gid(i-1,j-1,k) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k);
                mat[nelems] = Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i,j-1,k) < gidmax) {
                cols[nelems] = gid(i,j-1,k);
                mat[nelems] = Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i+1,j-1,k) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k);
                mat[nelems] = Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i-1,j,k) < gidmax) {
                cols[nelems] = gid(i-1,j,k);
                mat[nelems] = Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i+1,j,k) < gidmax) {
                cols[nelems] = gid(i+1,j,k);
                mat[nelems] = Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i-1,j+1,k) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k);
                mat[nelems] = Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i,j+1,k) < gidmax) {
                cols[nelems] = gid(i,j+1,k);
                mat[nelems] = Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2];
                ++nelems;
            }

            if                (gid(i+1,j+1,k) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k);
                mat[nelems] = Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
                ++nelems;
            }

            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

#ifdef AMREX_USE_GPU
template <typename HypreInt, typename AtomicInt>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndtslap_fill_ijmatrix_gpu (const int ps, const int i, const int j, const int k,
                                  const int offset, Box const& ndbx,
                                  Array4<AtomicInt const> const& gid,
                                  Array4<int const> const& lid,
                                  HypreInt* const ncols, HypreInt* const cols, Real* const mat,
                                  GpuArray<Real,3> const& s,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    if (lid(i,j,k) >= 0)
    {
        const Real h00 = dxinv[0]*dxinv[0];
        const Real h01 = dxinv[0]*dxinv[1];
        const Real h11 = dxinv[1]*dxinv[1];

        constexpr auto gidmax = std::numeric_limits<AtomicInt>::max();

        if (offset == 0) {
            cols[ps] = gid(i,j,k);
            mat[ps] = Real(-4./3.)*h00*s[0] + Real(-4./3.)*h11*s[2];
            int nc = 1;
            if (gid(i-1,j-1,k) < gidmax) { ++nc; }
            if (gid(i  ,j-1,k) < gidmax) { ++nc; }
            if (gid(i+1,j-1,k) < gidmax) { ++nc; }
            if (gid(i-1,j  ,k) < gidmax) { ++nc; }
            if (gid(i+1,j  ,k) < gidmax) { ++nc; }
            if (gid(i-1,j+1,k) < gidmax) { ++nc; }
            if (gid(i  ,j+1,k) < gidmax) { ++nc; }
            if (gid(i+1,j+1,k) < gidmax) { ++nc; }
            ncols[lid(i,j,k)] = nc;
        }
        else if (offset == 1 && gid(i-1,j-1,k) < gidmax) {
            cols[ps] =          gid(i-1,j-1,k);
            mat[ps] = Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
        }
        else if (offset == 2 && gid(i  ,j-1,k) < gidmax) {
            cols[ps] =          gid(i  ,j-1,k);
            mat[ps] = Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2];
        }
        else if (offset == 3 && gid(i+1,j-1,k) < gidmax) {
            cols[ps] =          gid(i+1,j-1,k);
            mat[ps] = Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
        }
        else if (offset == 4 && gid(i-1,j  ,k) < gidmax) {
            cols[ps] =          gid(i-1,j  ,k);
            mat[ps] = Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2];
        }
        else if (offset == 5 && gid(i+1,j  ,k) < gidmax) {
            cols[ps] =          gid(i+1,j  ,k);
            mat[ps] = Real(2./3.)*h00*s[0] - Real(1./3.)*h11*s[2];
        }
        else if (offset == 6 && gid(i-1,j+1,k) < gidmax) {
            cols[ps] =          gid(i-1,j+1,k);
            mat[ps] = Real(1./6.)*h00*s[0] - Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
        }
        else if (offset == 7 && gid(i  ,j+1,k) < gidmax) {
            cols[ps] =          gid(i  ,j+1,k);
            mat[ps] = Real(-1./3.)*h00*s[0] + Real(2./3.)*h11*s[2];
        }
        else if (offset == 8 && gid(i+1,j+1,k) < gidmax) {
            cols[ps] =          gid(i+1,j+1,k);
            mat[ps] = Real(1./6.)*h00*s[0] + Real(0.5)*h01*s[1] + Real(1./6.)*h11*s[2];
        }
    }
}
#endif

#endif

}

#endif
