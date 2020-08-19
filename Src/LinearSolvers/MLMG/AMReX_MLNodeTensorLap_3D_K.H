#ifndef AMREX_MLNODETENSORLAP_3D_K_H_
#define AMREX_MLNODETENSORLAP_3D_K_H_

namespace amrex {

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_x (Array4<Real const> const& crse,
                           int ic, int jc, int kc) noexcept
    {
        return (crse(ic,jc,kc)+crse(ic+1,jc,kc))*Real(0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_y (Array4<Real const> const& crse,
                           int ic, int jc, int kc) noexcept
    {
        return (crse(ic,jc,kc)+crse(ic,jc+1,kc))*Real(0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_line_z (Array4<Real const> const& crse,
                           int ic, int jc, int kc) noexcept
    {
        return (crse(ic,jc,kc)+crse(ic,jc,kc+1))*Real(0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_face_xy (Array4<Real const> const& crse,
                            int ic, int jc, int kc) noexcept
    {
        return (ts_interp_line_y(crse,ic  ,jc  ,kc) +
                ts_interp_line_y(crse,ic+1,jc  ,kc) +
                ts_interp_line_x(crse,ic  ,jc  ,kc) +
                ts_interp_line_x(crse,ic  ,jc+1,kc)) * Real(0.25);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_face_xz (Array4<Real const> const& crse,
                            int ic, int jc, int kc) noexcept
    {
        return (ts_interp_line_z(crse,ic  ,jc,kc  ) +
                ts_interp_line_z(crse,ic+1,jc,kc  ) +
                ts_interp_line_x(crse,ic  ,jc,kc  ) +
                ts_interp_line_x(crse,ic  ,jc,kc+1)) * Real(0.25);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ts_interp_face_yz (Array4<Real const> const& crse,
                            int ic, int jc, int kc) noexcept
    {
        return (ts_interp_line_z(crse,ic,jc  ,kc  ) +
                ts_interp_line_z(crse,ic,jc+1,kc  ) +
                ts_interp_line_y(crse,ic,jc  ,kc  ) +
                ts_interp_line_y(crse,ic,jc  ,kc+1)) * Real(0.25);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int i, int j, int k, Array4<Real> const& fine,
                          Array4<Real const> const& crse, Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd and j_is_odd and k_is_odd) {
            // Fine node at center of cell
            fine(i,j,k) += (ts_interp_face_yz(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_yz(crse,ic+1,jc  ,kc  ) +
                            ts_interp_face_xz(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_xz(crse,ic  ,jc+1,kc  ) +
                            ts_interp_face_xy(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_xy(crse,ic  ,jc  ,kc+1)) * Real(1./6.);
        } else if (j_is_odd and k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += ts_interp_face_yz(crse,ic,jc,kc);
        } else if (i_is_odd and k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += ts_interp_face_xz(crse,ic,jc,kc);
        } else if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += ts_interp_face_xy(crse,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += ts_interp_line_x(crse,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += ts_interp_line_y(crse,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += ts_interp_line_z(crse,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const& b, Array4<Real> const& y, Array4<Real const> const& x,
                      GpuArray<Real,6> const& s,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h02 = dxinv[0]*dxinv[2];
    const Real h11 = dxinv[1]*dxinv[1];
    const Real h12 = dxinv[1]*dxinv[2];
    const Real h22 = dxinv[2]*dxinv[2];
    amrex::LoopConcurrent(b, [=] (int i, int j, int k) noexcept
    {
        y(i,j,k) = (1./36.) *
            ( x(i-1,j-1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]))
            + x(i+1,j-1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]))
            + x(i-1,j+1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]))
            + x(i+1,j+1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]))
            + x(i-1,j-1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]))
            + x(i+1,j-1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]))
            + x(i-1,j+1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]))
            + x(i+1,j+1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]))
            + x(i-1,j-1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1])
            + x(i+1,j-1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1])
            + x(i-1,j+1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1])
            + x(i+1,j+1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1])
            + x(i-1,j  ,k-1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2])
            + x(i+1,j  ,k-1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2])
            + x(i-1,j  ,k+1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2])
            + x(i+1,j  ,k+1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2])
            + x(i  ,j-1,k-1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4])
            + x(i  ,j+1,k-1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4])
            + x(i  ,j-1,k+1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4])
            + x(i  ,j+1,k+1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4])
            + x(i-1,j  ,k  ) * (8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]))
            + x(i+1,j  ,k  ) * (8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]))
            + x(i  ,j-1,k  ) * (8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]))
            + x(i  ,j+1,k  ) * (8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]))
            + x(i  ,j  ,k-1) * (8.*(2.*h22*s[5]-h00*s[0]-h11*s[3]))
            + x(i  ,j  ,k+1) * (8.*(2.*h22*s[5]-h00*s[0]-h11*s[3]))
            + x(i  ,j  ,k  ) * (-32.)*(h00*s[0]+h11*s[3]+h22*s[5]) );
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_gauss_seidel (Box const& b, Array4<Real> const& sol,
                             Array4<Real const> const& rhs, Array4<int const> const& msk,
                             GpuArray<Real,6> const& s,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h02 = dxinv[0]*dxinv[2];
    const Real h11 = dxinv[1]*dxinv[1];
    const Real h12 = dxinv[1]*dxinv[2];
    const Real h22 = dxinv[2]*dxinv[2];
    amrex::LoopConcurrent(b, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (1./36.) * (-32.)*(h00*s[0]+h11*s[3]+h22*s[5]);
            Real Ax = (1./36.) *
                ( sol(i-1,j-1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]))
                + sol(i+1,j-1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]))
                + sol(i-1,j+1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]))
                + sol(i+1,j+1,k-1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]))
                + sol(i-1,j-1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]))
                + sol(i+1,j-1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]))
                + sol(i-1,j+1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]))
                + sol(i+1,j+1,k+1) * ((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]))
                + sol(i-1,j-1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1])
                + sol(i+1,j-1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1])
                + sol(i-1,j+1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1])
                + sol(i+1,j+1,k  ) * (4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1])
                + sol(i-1,j  ,k-1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2])
                + sol(i+1,j  ,k-1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2])
                + sol(i-1,j  ,k+1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2])
                + sol(i+1,j  ,k+1) * (4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2])
                + sol(i  ,j-1,k-1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4])
                + sol(i  ,j+1,k-1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4])
                + sol(i  ,j-1,k+1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4])
                + sol(i  ,j+1,k+1) * (4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4])
                + sol(i-1,j  ,k  ) * (8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]))
                + sol(i+1,j  ,k  ) * (8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]))
                + sol(i  ,j-1,k  ) * (8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]))
                + sol(i  ,j+1,k  ) * (8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]))
                + sol(i  ,j  ,k-1) * (8.*(2.*h22*s[5]-h00*s[0]-h11*s[3]))
                + sol(i  ,j  ,k+1) * (8.*(2.*h22*s[5]-h00*s[0]-h11*s[3])) )
                + sol(i  ,j  ,k  ) * s0;
            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_normalize (Box const& b, Array4<Real> const& phi,
                          Array4<int const> const& msk, GpuArray<Real,6> const& s,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h11 = dxinv[1]*dxinv[1];
    const Real h22 = dxinv[2]*dxinv[2];
    amrex::LoopConcurrent(b, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            Real s0 = (1./36.) * (-32.)*(h00*s[0]+h11*s[3]+h22*s[5]);
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
                              Array4<int const> const& msk, GpuArray<Real,6> const& s,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real h00 = dxinv[0]*dxinv[0];
    const Real h01 = dxinv[0]*dxinv[1];
    const Real h02 = dxinv[0]*dxinv[2];
    const Real h11 = dxinv[1]*dxinv[1];
    const Real h12 = dxinv[1]*dxinv[2];
    const Real h22 = dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(ndbx);
    const auto hi = amrex::ubound(ndbx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (nid(i,j,k) >= 0 and owner(i,j,k))
        {
            rows.push_back(nid(i,j,k));
            cols.push_back(nid(i,j,k));
            mat.push_back((1./36.) * (-32.)*(h00*s[0]+h11*s[3]+h22*s[5]));
            HypreNodeLap::Int nc = 1;

            if (               nid(i-1,j-1,k-1) >= 0) {
                cols.push_back(nid(i-1,j-1,k-1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j-1,k-1) >= 0) {
                cols.push_back(nid(i+1,j-1,k-1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j+1,k-1) >= 0) {
                cols.push_back(nid(i-1,j+1,k-1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j+1,k-1) >= 0) {
                cols.push_back(nid(i+1,j+1,k-1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j-1,k+1) >= 0) {
                cols.push_back(nid(i-1,j-1,k+1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]-h02*s[2]-h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j-1,k+1) >= 0) {
                cols.push_back(nid(i+1,j-1,k+1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]+h02*s[2]-h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j+1,k+1) >= 0) {
                cols.push_back(nid(i-1,j+1,k+1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*(-h01*s[1]-h02*s[2]+h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j+1,k+1) >= 0) {
                cols.push_back(nid(i+1,j+1,k+1));
                Real t = (1./36.)*((h00*s[0]+h11*s[3]+h22*s[5])+3.*( h01*s[1]+h02*s[2]+h12*s[4]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j-1,k  ) >= 0) {
                cols.push_back(nid(i-1,j-1,k  ));
                Real t = (1./36.)*(4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j-1,k  ) >= 0) {
                cols.push_back(nid(i+1,j-1,k  ));
                Real t = (1./36.)*(4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j+1,k  ) >= 0) {
                cols.push_back(nid(i-1,j+1,k  ));
                Real t = (1./36.)*(4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])-12.*h01*s[1]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j+1,k  ) >= 0) {
                cols.push_back(nid(i+1,j+1,k  ));
                Real t = (1./36.)*(4.*(h00*s[0]+h11*s[3]-0.5*h22*s[5])+12.*h01*s[1]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j  ,k-1) >= 0) {
                cols.push_back(nid(i-1,j  ,k-1));
                Real t = (1./36.)*(4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j  ,k-1) >= 0) {
                cols.push_back(nid(i+1,j  ,k-1));
                Real t = (1./36.)*(4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j  ,k+1) >= 0) {
                cols.push_back(nid(i-1,j  ,k+1));
                Real t = (1./36.)*(4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])-12.*h02*s[2]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j  ,k+1) >= 0) {
                cols.push_back(nid(i+1,j  ,k+1));
                Real t = (1./36.)*(4.*(h00*s[0]+h22*s[5]-0.5*h11*s[3])+12.*h02*s[2]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j-1,k-1) >= 0) {
                cols.push_back(nid(i  ,j-1,k-1));
                Real t = (1./36.)*(4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j+1,k-1) >= 0) {
                cols.push_back(nid(i  ,j+1,k-1));
                Real t = (1./36.)*(4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j-1,k+1) >= 0) {
                cols.push_back(nid(i  ,j-1,k+1));
                Real t = (1./36.)*(4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])-12.*h12*s[4]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j+1,k+1) >= 0) {
                cols.push_back(nid(i  ,j+1,k+1));
                Real t = (1./36.)*(4.*(h11*s[3]+h22*s[5]-0.5*h00*s[0])+12.*h12*s[4]);
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i-1,j  ,k  ) >= 0) {
                cols.push_back(nid(i-1,j  ,k  ));
                Real t = (1./36.)*(8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i+1,j  ,k  ) >= 0) {
                cols.push_back(nid(i+1,j  ,k  ));
                Real t = (1./36.)*(8.*(2.*h00*s[0]-h11*s[3]-h22*s[5]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j-1,k  ) >= 0) {
                cols.push_back(nid(i  ,j-1,k  ));
                Real t = (1./36.)*(8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j+1,k  ) >= 0) {
                cols.push_back(nid(i  ,j+1,k  ));
                Real t = (1./36.)*(8.*(2.*h11*s[3]-h00*s[0]-h22*s[5]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j  ,k-1) >= 0) {
                cols.push_back(nid(i  ,j  ,k-1));
                Real t = (1./36.)*(8.*(2.*h22*s[5]-h00*s[0]-h11*s[3]));
                mat.push_back(t);
                ++nc;
            }

            if (               nid(i  ,j  ,k+1) >= 0) {
                cols.push_back(nid(i  ,j  ,k+1));
                Real t = (1./36.)*(8.*(2.*h22*s[5]-h00*s[0]-h11*s[3]));
                mat.push_back(t);
                ++nc;
            }

            ncols.push_back(nc);
        }
    }}}
}

#endif

}

#endif
