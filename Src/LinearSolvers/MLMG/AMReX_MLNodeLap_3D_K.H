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

//
// coeffs
//

//
// bc
//

//
// operator
//

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
