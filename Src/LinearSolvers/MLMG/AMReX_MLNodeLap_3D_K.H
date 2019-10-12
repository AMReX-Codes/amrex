#ifndef AMREX_MLNODELAP_3D_K_H_
#define AMREX_MLNODELAP_3D_K_H_

namespace amrex {

//
// masks
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int i, int j, int k, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    int s = cmsk(i-1,j-1,k-1) + cmsk(i  ,j-1,k-1)
        +   cmsk(i-1,j  ,k-1) + cmsk(i  ,j  ,k-1)
        +   cmsk(i-1,j-1,k  ) + cmsk(i  ,j-1,k  )
        +   cmsk(i-1,j  ,k  ) + cmsk(i  ,j  ,k  );
    if (s == 8*crse_cell) {
        nmsk(i,j,k) = crse_node;
    }
    else if (s == 8*fine_cell) {
        nmsk(i,j,k) = fine_node;
    } else {
        nmsk(i,j,k) = crse_fine_node;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 Array<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 Array<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,k) = (omsk(i-1,j-1,k-1) == 1 or omsk(i,j-1,k-1) == 1 or
                       omsk(i-1,j  ,k-1) == 1 or omsk(i,j  ,k-1) == 1 or
                       omsk(i-1,j-1,k  ) == 1 or omsk(i,j-1,k  ) == 1 or
                       omsk(i-1,j  ,k  ) == 1 or omsk(i,j  ,k  ) == 1);
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if (bclo[0] == LinOpBCType::Dirichlet and lo.x == domlo.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) = 1;
        }}
    }

    if (bchi[0] == LinOpBCType::Dirichlet and hi.x == domhi.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) = 1;
        }}
    }

    if (bclo[1] == LinOpBCType::Dirichlet and lo.y == domlo.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) = 1;
        }}
    }

    if (bchi[1] == LinOpBCType::Dirichlet and hi.y == domhi.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) = 1;
        }}
    }

    if (bclo[2] == LinOpBCType::Dirichlet and lo.z == domlo.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) = 1;
        }}
    }

    if (bchi[2] == LinOpBCType::Dirichlet and hi.z == domhi.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) = 1;
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           Array<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           Array<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,k) = static_cast<Real>(omsk(i,j,k));
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
        and lo.x == domlo.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) *= 0.5;
        }}
    }

    if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
        and hi.x == domhi.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) *= 0.5;
        }}
    }

    if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
        and lo.y == domlo.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) *= 0.5;
        }}
    }

    if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
        and hi.y == domhi.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) *= 0.5;
        }}
    }

    if ((bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)
        and lo.z == domlo.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) *= 0.5;
        }}
    }

    if ((bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)
        and hi.z == domhi.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) *= 0.5;
        }}
    }
}

//
// coeffs
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i  ,2*j,2*k  )+fine(2*i  ,2*j+1,2*k  )+
              fine(2*i  ,2*j,2*k+1)+fine(2*i  ,2*j+1,2*k+1);
    Real cr = fine(2*i+1,2*j,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i+1,2*j,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_y (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1);
    Real cr = fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_z (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  );
    Real cr = fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1)+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

//
// bc
//

//
// operator
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_ha (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<Real const> const& sy,
                       Array4<Real const> const& sz, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1./36.)*dxinv[0]*dxinv[0];
    Real facy = (1./36.)*dxinv[1]*dxinv[1];
    Real facz = (1./36.)*dxinv[2]*dxinv[2];
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            y(i,j,k) = 0.0;
        } else {
            y(i,j,k) = x(i,j,k)*(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                             +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                       +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                             +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                       +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                             +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
            y(i,j,k) += x(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                       +facy*sy(i-1,j-1,k-1)
                                       +facz*sz(i-1,j-1,k-1))
                      + x(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                       +facy*sy(i  ,j-1,k-1)
                                       +facz*sz(i  ,j-1,k-1))
                      + x(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                       +facy*sy(i-1,j  ,k-1)
                                       +facz*sz(i-1,j  ,k-1))
                      + x(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                       +facy*sy(i  ,j  ,k-1)
                                       +facz*sz(i  ,j  ,k-1))
                      + x(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                       +facy*sy(i-1,j-1,k  )
                                       +facz*sz(i-1,j-1,k  ))
                      + x(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                       +facy*sy(i  ,j-1,k  )
                                       +facz*sz(i  ,j-1,k  ))
                      + x(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                       +facy*sy(i-1,j  ,k  )
                                       +facz*sz(i-1,j  ,k  ))
                      + x(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                       +facy*sy(i  ,j  ,k  )
                                       +facz*sz(i  ,j  ,k  ));
            y(i,j,k) += x(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                      + x(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                        +2.0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                      + x(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                        +2.0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                      + x(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                        +2.0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                        +2.0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                      + x(i-1,j  ,k-1)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                      + x(i+1,j  ,k-1)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                        +2.0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                      + x(i-1,j  ,k+1)*( 2.0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                      + x(i+1,j  ,k+1)*( 2.0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                            -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                        +2.0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                      + x(i-1,j-1,k  )*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                      + x(i+1,j-1,k  )*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                        +2.0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                      + x(i-1,j+1,k  )*( 2.0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                      + x(i+1,j+1,k  )*( 2.0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                        +2.0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                            -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)));
                y(i,j,k) += 2.0*x(i-1,j,k)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                                -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                                -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                          + 2.0*x(i+1,j,k)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                                -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                                -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                          + 2.0*x(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                            +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                                -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                          + 2.0*x(i,j+1,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                            +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                                -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                          + 2.0*x(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                                -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                            +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                          + 2.0*x(i,j,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                                -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                            +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_aa (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + 2.0*facy + 2.0*facz;
    Real f2xmy2z = 2.0*facx - facy + 2.0*facz;
    Real f2x2ymz = 2.0*facx + 2.0*facy - facz;
    Real f4xm2ym2z = 4.0*facx - 2.0*facy - 2.0*facz;
    Real fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
    Real fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            y(i,j,k) = 0.0;
        } else {
            y(i,j,k) = x(i,j,k)*(-4.0)*fxyz*
                (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
                + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                      + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                      + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                      + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                      + x(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                      + x(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                      + x(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                      + x(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                         + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                         + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                         + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                         + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                         + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                         + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                         + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                         + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                         + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                           + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                           + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                           + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sx,
                           Array4<Real const> const& sy, Array4<Real const> const& sz,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) = x(i,j,k)/((-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                              +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                        +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                              +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                              +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  ))));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sig,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) = x(i,j,k) /
                ((-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                             +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sx,
                        Array4<Real const> const& sy, Array4<Real const> const& sz,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = -4.0 * (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = -4.0 * (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = -4.0 * (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                        +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                  +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                        +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                  +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                        +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sig,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real fxyz = -4.0 * (1.0/36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                        +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sx,
                              Array4<Real const> const& sy, Array4<Real const> const& sz,
                              Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                   +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                             +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                   +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                             +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                   +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
            Real Ax = sol(i,j,k)*s0
                     + sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                        +facy*sy(i-1,j-1,k-1)
                                        +facz*sz(i-1,j-1,k-1))
                     + sol(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                        +facy*sy(i  ,j-1,k-1)
                                        +facz*sz(i  ,j-1,k-1))
                     + sol(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                        +facy*sy(i-1,j  ,k-1)
                                        +facz*sz(i-1,j  ,k-1))
                     + sol(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                        +facy*sy(i  ,j  ,k-1)
                                        +facz*sz(i  ,j  ,k-1))
                     + sol(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                        +facy*sy(i-1,j-1,k  )
                                        +facz*sz(i-1,j-1,k  ))
                     + sol(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                        +facy*sy(i  ,j-1,k  )
                                        +facz*sz(i  ,j-1,k  ))
                     + sol(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                        +facy*sy(i-1,j  ,k  )
                                        +facz*sz(i-1,j  ,k  ))
                     + sol(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                        +facy*sy(i  ,j  ,k  )
                                        +facz*sz(i  ,j  ,k  ))
                     +sol(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                     +sol(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                        +2.0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                     +sol(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                        +2.0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                     +sol(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                        +2.0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                        +2.0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                     +sol(i-1,j  ,k-1)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                     +sol(i+1,j  ,k-1)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                        +2.0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                     +sol(i-1,j  ,k+1)*( 2.0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                     +sol(i+1,j  ,k+1)*( 2.0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                            -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                        +2.0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                     +sol(i-1,j-1,k  )*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                     +sol(i+1,j-1,k  )*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                        +2.0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                     +sol(i-1,j+1,k  )*( 2.0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                     +sol(i+1,j+1,k  )*( 2.0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                        +2.0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                            -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                     + 2.0*sol(i-1,j,k)*(2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                     + 2.0*sol(i+1,j,k)*(2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                     + 2.0*sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                     + 2.0*sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                     + 2.0*sol(i,j,k-1)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                     + 2.0*sol(i,j,k+1)*(   -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));

                sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + 2.0*facy + 2.0*facz;
    Real f2xmy2z = 2.0*facx - facy + 2.0*facz;
    Real f2x2ymz = 2.0*facx + 2.0*facy - facz;
    Real f4xm2ym2z = 4.0*facx - 2.0*facy - 2.0*facz;
    Real fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
    Real fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                  +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
            Real Ax = sol(i,j,k)*s0
                + fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                         + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                         + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                         + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                         + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                         + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                         + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                         + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                         + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                         + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                           + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                           + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                           + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

//
// restriction
//

//
// interpolation
//

//
// rhs & u
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.25*dxinv[0];
    Real facy = 0.25*dxinv[1];
    Real facz = 0.25*dxinv[2];

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = facx*(-vel(i-1,j-1,k-1,0)+vel(i,j-1,k-1,0)
                           -vel(i-1,j  ,k-1,0)+vel(i,j  ,k-1,0)
                           -vel(i-1,j-1,k  ,0)+vel(i,j-1,k  ,0)
                           -vel(i-1,j  ,k  ,0)+vel(i,j  ,k  ,0))
                   + facy*(-vel(i-1,j-1,k-1,1)-vel(i,j-1,k-1,1)
                           +vel(i-1,j  ,k-1,1)+vel(i,j  ,k-1,1)
                           -vel(i-1,j-1,k  ,1)-vel(i,j-1,k  ,1)
                           +vel(i-1,j  ,k  ,1)+vel(i,j  ,k  ,1))
                   + facz*(-vel(i-1,j-1,k-1,2)-vel(i,j-1,k-1,2)
                           -vel(i-1,j  ,k-1,2)-vel(i,j  ,k-1,2)
                           +vel(i-1,j-1,k  ,2)+vel(i,j-1,k  ,2)
                           +vel(i-1,j  ,k  ,2)+vel(i,j  ,k  ,2));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = 0.125 * (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) +
                              rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, Real dxi, Real dyi, Real dzi) noexcept
{
    Real facx = 0.25*dxi;
    Real facy = 0.25*dyi;
    Real facz = 0.25*dzi;
    u(i,j,k,0) -= sig(i,j,k)*facx
        * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,1) -= sig(i,j,k)*facy
        * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,2) -= sig(i,j,k)*facz
        * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
           +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
}

//
// residual
//

//
// sync residual
//

}

#endif
