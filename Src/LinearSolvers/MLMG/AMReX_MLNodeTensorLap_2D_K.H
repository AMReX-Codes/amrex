#ifndef AMREX_MLNODETENSORLAP_2D_K_H_
#define AMREX_MLNODETENSORLAP_2D_K_H_

namespace amrex {

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_x (Array4<Real const> const& crse, int i, int j, int ic, int jc) noexcept
    {
        return (crse(ic,jc,0)+crse(ic+1,jc,0))*0.5_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_y (Array4<Real const> const& crse, int i, int j, int ic, int jc) noexcept
    {
        return (crse(ic,jc,0)+crse(ic,jc+1,0))*0.5_rt;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_face_xy (Array4<Real const> const& crse, int i, int j, int ic, int jc) noexcept
    {
        return (ts_interp_line_y(crse,i-1,j  ,ic  ,jc  ) +
                ts_interp_line_y(crse,i+1,j  ,ic+1,jc  ) +
                ts_interp_line_x(crse,i  ,j-1,ic  ,jc  ) +
                ts_interp_line_x(crse,i  ,j+1,ic  ,jc+1)) * 0.25_rt;
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
            fine(i,j,0) += ts_interp_face_xy(crse,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += ts_interp_line_x(crse,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += ts_interp_line_y(crse,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const& b, Array4<Real> const& y, Array4<Real const> const& x,
                      Array4<int const> const& msk, GpuArray<Real,3> const& s,
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

}

#endif
