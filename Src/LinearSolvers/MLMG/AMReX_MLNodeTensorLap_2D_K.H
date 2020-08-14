#ifndef AMREX_MLNODETENSORLAP_2D_K_H_
#define AMREX_MLNODETENSORLAP_2D_K_H_

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
        if (i_is_odd and j_is_odd) {
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
        y(i,j,k) = x(i-1,j-1,0) * ((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2])
            +      x(i-1,j  ,0) * ((2./3.)*h00*s[0] - (1./3.)*h11*s[2])
            +      x(i-1,j+1,0) * ((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2])
            +      x(i  ,j-1,0) * ((-1./3.)*h00*s[0] + (2./3.)*h11*s[2])
            +      x(i  ,j  ,0) * ((-4./3.)*h00*s[0] + (-4./3.)*h11*s[2])
            +      x(i  ,j+1,0) * ((-1./3.)*h00*s[0] + (2./3.)*h11*s[2])
            +      x(i+1,j-1,0) * ((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2])
            +      x(i+1,j  ,0) * ((2./3.)*h00*s[0] - (1./3.)*h11*s[2])
            +      x(i+1,j+1,0) * ((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2]);
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
            Real s0 = ((-4./3.)*h00*s[0] + (-4./3.)*h11*s[2]);
            Real Ax = sol(i-1,j-1,0) * ((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2])
                +     sol(i-1,j  ,0) * ((2./3.)*h00*s[0] - (1./3.)*h11*s[2])
                +     sol(i-1,j+1,0) * ((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2])
                +     sol(i  ,j-1,0) * ((-1./3.)*h00*s[0] + (2./3.)*h11*s[2])
                +     sol(i  ,j  ,0) * s0
                +     sol(i  ,j+1,0) * ((-1./3.)*h00*s[0] + (2./3.)*h11*s[2])
                +     sol(i+1,j-1,0) * ((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2])
                +     sol(i+1,j  ,0) * ((2./3.)*h00*s[0] - (1./3.)*h11*s[2])
                +     sol(i+1,j+1,0) * ((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2]);
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
            Real s0 = ((-4./3.)*h00*s[0] + (-4./3.)*h11*s[2]);
            phi(i,j,k) /= s0;
        }
    });
}

#ifdef AMREX_USE_HYPRE

AMREX_FORCE_INLINE
void mlndtslap_fill_ijmatrix (Box const& ndbx, Array4<HypreNodeLap::Int const> const& nid,
                              Array4<int const> const& owner,
                              Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                              Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat,
                              Array4<int const> const& msk, GpuArray<Real,3> const& s,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h11 = dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(ndbx);
    const auto hi = amrex::ubound(ndbx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (nid(i,j,k) >= 0 and owner(i,j,k))
        {
            rows.push_back(nid(i,j,k));
            cols.push_back(nid(i,j,k));
            mat.push_back((-4./3.)*h00*s[0] + (-4./3.)*h11*s[2]);
            HypreNodeLap::Int nc = 1;

            if                (nid(i-1,j-1,k) >= 0) {
                cols.push_back(nid(i-1,j-1,k));
                mat.push_back((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i-1,j,k) >= 0) {
                cols.push_back(nid(i-1,j,k));
                mat.push_back((2./3.)*h00*s[0] - (1./3.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i-1,j+1,k) >= 0) {
                cols.push_back(nid(i-1,j+1,k));
                mat.push_back((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i,j-1,k) >= 0) {
                cols.push_back(nid(i,j-1,k));
                mat.push_back((-1./3.)*h00*s[0] + (2./3.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i,j+1,k) >= 0) {
                cols.push_back(nid(i,j+1,k));
                mat.push_back((-1./3.)*h00*s[0] + (2./3.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i+1,j-1,k) >= 0) {
                cols.push_back(nid(i+1,j-1,k));
                mat.push_back((1./6.)*h00*s[0] - 0.5*h01*s[1] + (1./6.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i+1,j,k) >= 0) {
                cols.push_back(nid(i+1,j,k));
                mat.push_back((2./3.)*h00*s[0] - (1./3.)*h11*s[2]);
                ++nc;
            }

            if                (nid(i+1,j+1,k) >= 0) {
                cols.push_back(nid(i+1,j+1,k));
                mat.push_back((1./6.)*h00*s[0] + 0.5*h01*s[1] + (1./6.)*h11*s[2]);
                ++nc;
            }

            ncols.push_back(nc);
        }
    }}}
}

#endif

}

#endif
