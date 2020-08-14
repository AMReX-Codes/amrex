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
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if (!dmsk(i,j,k)) {
            dmsk(i,j,k) = (omsk(i-1,j-1,k-1) == 1 or omsk(i,j-1,k-1) == 1 or
                           omsk(i-1,j  ,k-1) == 1 or omsk(i,j  ,k-1) == 1 or
                           omsk(i-1,j-1,k  ) == 1 or omsk(i,j-1,k  ) == 1 or
                           omsk(i-1,j  ,k  ) == 1 or omsk(i,j  ,k  ) == 1);
        }
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
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int i, int j, int k, Array4<Real> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,k-1) == fine_flag and
        msk(i  ,j-1,k-1) == fine_flag and
        msk(i-1,j  ,k-1) == fine_flag and
        msk(i  ,j  ,k-1) == fine_flag and
        msk(i-1,j-1,k  ) == fine_flag and
        msk(i  ,j-1,k  ) == fine_flag and
        msk(i-1,j  ,k  ) == fine_flag and
        msk(i  ,j  ,k  ) == fine_flag)
    {
        phi(i,j,k) = 0.0;
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine, int idir) noexcept
{
    if (idir == 2) {
        Real cl = fine(2*i  ,2*j,k) + fine(2*i  ,2*j+1,k);
        Real cr = fine(2*i+1,2*j,k) + fine(2*i+1,2*j+1,k);
        crse(i,j,k) = cl*cr/(cl+cr);    
    } else if (idir == 1) {
        Real cl = fine(2*i  ,j,2*k) + fine(2*i  ,j,2*k+1);
        Real cr = fine(2*i+1,j,2*k) + fine(2*i+1,j,2*k+1);
        crse(i,j,k) = cl*cr/(cl+cr);
    } else {
        Real cl = fine(i,2*j  ,2*k) + fine(i,2*j  ,2*k+1);
        Real cr = fine(i,2*j+1,2*k) + fine(i,2*j+1,2*k+1);
        crse(i,j,k) = cl*cr/(cl+cr);
    }
}
//
// bc
//

template <typename T>
inline void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                             GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                             GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    Box gdomain = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (not bflo[idim]) gdomain.growLo(idim,1);
        if (not bfhi[idim]) gdomain.growHi(idim,1);
    }

    if (gdomain.strictly_contains(vbx)) return;

    int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Box const& sbox = amrex::grow(vbx,1);
    AMREX_HOST_DEVICE_FOR_3D(sbox, i, j, k,
    {
        if (not gdomain.contains(IntVect(i,j,k))) {
            // xlo & ylo & zlo
            if (i == dlo.x-1 and j == dlo.y-1 and k == dlo.z-1 and (bflo[0] or bflo[1] or bflo[2]))
            {
                if (bflo[0] and bflo[1] and bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k+1+offset);
                }
                else if (bflo[0] and bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0] and bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bflo[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & ylo & zlo
            else if (i == dhi.x+1 and j == dlo.y-1 and k == dlo.z-1 and (bfhi[0] or bflo[1] or bflo[2]))
            {
                if (bfhi[0] and bflo[1] and bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k+1+offset);
                }
                else if (bfhi[0] and bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0] and bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bflo[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & yhi & zlo
            else if (i == dlo.x-1 and j == dhi.y+1 and k == dlo.z-1 and (bflo[0] or bfhi[1] or bflo[2]))
            {
                if (bflo[0] and bfhi[1] and bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k+1+offset);
                }
                else if (bflo[0] and bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0] and bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bfhi[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & yhi & zlo
            else if (i == dhi.x+1 and j == dhi.y+1 and k == dlo.z-1 and (bfhi[0] or bfhi[1] or bflo[2]))
            {
                if (bfhi[0] and bfhi[1] and bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k+1+offset);
                }
                else if (bfhi[0] and bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0] and bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bfhi[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & ylo & zhi
            else if (i == dlo.x-1 and j == dlo.y-1 and k == dhi.z+1 and (bflo[0] or bflo[1] or bfhi[2]))
            {
                if (bflo[0] and bflo[1] and bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k-1-offset);
                }
                else if (bflo[0] and bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0] and bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bflo[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & ylo & zhi
            else if (i == dhi.x+1 and j == dlo.y-1 and k == dhi.z+1 and (bfhi[0] or bflo[1] or bfhi[2]))
            {
                if (bfhi[0] and bflo[1] and bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k-1-offset);
                }
                else if (bfhi[0] and bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0] and bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bflo[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xlo & yhi & zhi
            else if (i == dlo.x-1 and j == dhi.y+1 and k == dhi.z+1 and (bflo[0] or bfhi[1] or bfhi[2]))
            {
                if (bflo[0] and bfhi[1] and bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k-1-offset);
                }
                else if (bflo[0] and bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0] and bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bfhi[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & yhi & zhi
            else if (i == dhi.x+1 and j == dhi.y+1 and k == dhi.z+1 and (bfhi[0] or bfhi[1] or bfhi[2]))
            {
                if (bfhi[0] and bfhi[1] and bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k-1-offset);
                }
                else if (bfhi[0] and bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0] and bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bfhi[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xlo & ylo
            else if (i == dlo.x-1 and j == dlo.y-1 and (bflo[0] or bflo[1]))
            {
                if (bflo[0] and bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xhi & ylo
            else if (i == dhi.x+1 and j == dlo.y-1 and (bfhi[0] or bflo[1]))
            {
                if (bfhi[0] and bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xlo & yhi
            else if (i == dlo.x-1 and j == dhi.y+1 and (bflo[0] or bfhi[1]))
            {
                if (bflo[0] and bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            // xhi & yhi
            else if (i == dhi.x+1 and j == dhi.y+1 and (bfhi[0] or bfhi[1]))
            {
                if (bfhi[0] and bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            // xlo & zlo
            else if (i == dlo.x-1 and k == dlo.z-1 and (bflo[0] or bflo[2]))
            {
                if (bflo[0] and bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & zlo
            else if (i == dhi.x+1 and k == dlo.z-1 and (bfhi[0] or bflo[2]))
            {
                if (bfhi[0] and bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & zhi
            else if (i == dlo.x-1 and k == dhi.z+1 and (bflo[0] or bfhi[2]))
            {
                if (bflo[0] and bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & zhi
            else if (i == dhi.x+1 and k == dhi.z+1 and (bfhi[0] or bfhi[2]))
            {
                if (bfhi[0] and bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // ylo & zlo
            else if (j == dlo.y-1 and k == dlo.z-1 and (bflo[1] or bflo[2]))
            {
                if (bflo[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // yhi & zlo
            else if (j == dhi.y+1 and k == dlo.z-1 and (bfhi[1] or bflo[2]))
            {
                if (bfhi[1] and bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // ylo & zhi
            else if (j == dlo.y-1 and k == dhi.z+1 and (bflo[1] or bfhi[2]))
            {
                if (bflo[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // yhi & zhi
            else if (j == dhi.y+1 and k == dhi.z+1 and (bfhi[1] or bfhi[2]))
            {
                if (bfhi[1] and bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            else if (i == dlo.x-1 and bflo[0])
            {
                a(i,j,k) = a(i+1+offset, j, k);
            }
            else if (i == dhi.x+1 and bfhi[0])
            {
                a(i,j,k) = a(i-1-offset, j, k);
            }
            else if (j == dlo.y-1 and bflo[1])
            {
                a(i,j,k) = a(i, j+1+offset, k);
            }
            else if (j == dhi.y+1 and bfhi[1])
            {
                a(i,j,k) = a(i, j-1-offset, k);
            }
            else if (k == dlo.z-1 and bflo[2])
            {
                a(i,j,k) = a(i, j, k+1+offset);
            }
            else if (k == dhi.z+1 and bfhi[2])
            {
                a(i,j,k) = a(i, j, k-1-offset);
            }
        }
    });
}

//
// operator
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_ha (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<Real const> const& sy,
                       Array4<Real const> const& sz, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        Real facx = (1./36.)*dxinv[0]*dxinv[0];
        Real facy = (1./36.)*dxinv[1]*dxinv[1];
        Real facz = (1./36.)*dxinv[2]*dxinv[2];
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
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_aa (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tridiagonal_solve (Array1D<Real,0,31>& a_ls, Array1D<Real,0,31>& b_ls, Array1D<Real,0,31>& c_ls,
                        Array1D<Real,0,31>& r_ls, Array1D<Real,0,31>& u_ls, Array1D<Real,0,31>& gam,
                        int ilen ) noexcept
{
    Real bet = b_ls(0);
    u_ls(0) = r_ls(0) / bet;

    for (int i = 1; i <= ilen - 1; i++) {
        gam(i) = c_ls(i-1) / bet;
        bet = b_ls(i) - a_ls(i)*gam(i);
        if (bet == 0) amrex::Abort(">>>TRIDIAG FAILED");
        u_ls(i) = (r_ls(i)-a_ls(i)*u_ls(i-1)) / bet;
    }
    for (int i = ilen-2; i >= 0; i--) {
        u_ls(i) = u_ls(i) - gam(i+1)*u_ls(i+1);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_with_line_solve_aa (Box const& bx, Array4<Real> const& sol,
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
 
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    int idir = -1;
    int ilen = 33;

    if ( (dxinv[0] <= dxinv[2]) and (dxinv[1] <= dxinv[2]) ) {
        idir = 2;
        ilen = hi.z - lo.z + 1;
    }
    if ( (dxinv[0] <= dxinv[1]) and (dxinv[2] <= dxinv[1]) ) {
        idir = 1;
        ilen = hi.y - lo.y + 1;
    }
    if ( (dxinv[1] <= dxinv[0]) and (dxinv[2] <= dxinv[0]) ) {
        idir = 0;
        ilen = hi.x - lo.x + 1;
    }

    if (ilen > 32) amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is hard-wired to be no longer than 32");

    Array1D<Real,0,31> a_ls,b_ls,c_ls,u_ls,r_ls,gam;


    if ( idir == 2 ) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                for (int k = lo.z; k <= hi.z; ++k)
                {
                    if (msk(i,j,k)) 
                    {
                        a_ls(k-lo.z) = 0.;
                        b_ls(k-lo.z) = 1.;
                        c_ls(k-lo.z) = 0.;
                        u_ls(k-lo.z) = 0.;
                        r_ls(k-lo.z) = 0.;
                    } 
                    else 
                    {
                        Real s0 = (-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)   
                                             + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        Real Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
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
                                       + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)));

                        a_ls(k-lo.z) = fm2xm2y4z*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1));
                        b_ls(k-lo.z) = s0;
                        c_ls(k-lo.z) = fm2xm2y4z*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
                        u_ls(k-lo.z) = 0.;  
                        r_ls(k-lo.z) = rhs(i,j,k) - Ax;
                    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int k = lo.z; k <= hi.z; ++k)
                {  
                    sol(i,j,k) = u_ls(k-lo.z);
                }
            }
        }
    }
    else if (idir == 1)
    {
        for (int k = lo.z; k <= hi.z; ++k)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                for (int j = lo.y; j <= hi.y; ++j)
                {
                    if (msk(i,j,k)) {
                        a_ls(j-lo.y) = 0.;
                        b_ls(j-lo.y) = 1.;
                        c_ls(j-lo.y) = 0.;
                        u_ls(j-lo.y) = 0.;
                        r_ls(j-lo.y) = 0.;
                    } 
                    else
                    {
                        Real s0 = (-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                             + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        Real Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
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
                            + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                                       + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

                        a_ls(j-lo.y) = fm2x4ym2z*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k));
                        b_ls(j-lo.y) = s0;
                        c_ls(j-lo.y) = fm2x4ym2z*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k));
                        u_ls(j-lo.y) = 0.;
                        r_ls(j-lo.y) = rhs(i,j,k) - Ax;

                    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int j = lo.y; j <= hi.y; ++j)
                {
                    sol(i,j,k) = u_ls(j-lo.y);
                }
            }
        }
    } 
    else if (idir == 0) 
    {
        for (int j = lo.y; j <= hi.y; ++j)
        {
            for (int k = lo.z; k <= hi.z; ++k)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if (msk(i,j,k)) 
                    {
                        a_ls(i-lo.x) = 0.;
                        b_ls(i-lo.x) = 1.;
                        c_ls(i-lo.x) = 0.;
                        u_ls(i-lo.x) = 0.;
                        r_ls(i-lo.x) = 0.;
                    }
                    else
                    {
                        Real s0 = (-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                             + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        Real Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
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
                            + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                                     + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                            + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                                     + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

                        a_ls(i-lo.x) = f4xm2ym2z*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k));
                        b_ls(i-lo.x) = s0;
                        c_ls(i-lo.x) = f4xm2ym2z*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k));
                        u_ls(i-lo.x) = 0.;
                        r_ls(i-lo.x) = rhs(i,j,k) - Ax;
              	    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int i = lo.x; i <= hi.x; ++i)
                {
                    sol(i,j,k) = u_ls(i-lo.x);
                }
            }
        }
    } 
    else
    {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is wrong direction.");	
    } 
}

//
// restriction
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {
        crse(i,j,k) = (1./64.)*(fine(ii-1,jj-1,kk-1)+fine(ii+1,jj-1,kk-1)
                               +fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1)
                               +fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1)
                               +fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1))
                    + (1./32.)*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1)
                               +fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)
                               +fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1)
                               +fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)
                               +fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  )
                               +fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ))
                    + (1./16.)*(fine(ii-1,jj,kk)+fine(ii+1,jj,kk)
                               +fine(ii,jj-1,kk)+fine(ii,jj+1,kk)
                               +fine(ii,jj,kk-1)+fine(ii,jj,kk+1))
                      + (1./8.)*fine(ii,jj,kk);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int i, int j, int k, Array4<Real> const& crse,
                               Array4<Real const> const& fine, Array4<int const> const& msk, int idir) noexcept
{
    if (idir == 2) 
    {
        int ii = i*2;
        int jj = j*2;
        int kk = k;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = 0.0;
        } else { // use 2-D formula 
            crse(i,j,k) = (1./16.)*(fine(ii-1,jj-1,kk) + 2.*fine(ii  ,jj-1,kk) +    fine(ii+1,jj-1,kk)
                               + 2.*fine(ii-1,jj  ,kk) + 4.*fine(ii  ,jj  ,kk) + 2.*fine(ii+1,jj  ,kk)
                                  + fine(ii-1,jj+1,kk) + 2.*fine(ii  ,jj+1,kk) +    fine(ii+1,jj+1,kk));
        }
    } 
    else if (idir == 1) 
    {
        int ii = i*2;
        int jj = j;
        int kk = k*2;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = 0.0;
        } else { // use 2-D formula 
            crse(i,j,k) = (1./16.)*(fine(ii-1,jj,kk-1) + 2.*fine(ii  ,jj,kk-1) +    fine(ii+1,jj,kk-1)
                               + 2.*fine(ii-1,jj  ,kk) + 4.*fine(ii  ,jj,kk  ) + 2.*fine(ii+1,jj,kk  )
                                  + fine(ii-1,jj,kk+1) + 2.*fine(ii  ,jj,kk+1) +    fine(ii+1,jj,kk+1));
        }
    }
    else if (idir == 0)
    {
        int ii = i;
        int jj = j*2;
        int kk = k*2;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = 0.0;
        } else { // use 2-D formula 
            crse(i,j,k) = (1./16.)*(fine(ii,jj-1,kk-1) + 2.*fine(ii  ,jj,kk-1) +    fine(ii,jj+1,kk-1)
                               + 2.*fine(ii,jj-1  ,kk) + 4.*fine(ii  ,jj,kk  ) + 2.*fine(ii,jj+1,kk  )
                                  + fine(ii,jj-1,kk+1) + 2.*fine(ii  ,jj,kk+1) +    fine(ii,jj+1,kk+1));
        }
    }
    else 
    {
        amrex::Abort("mlndlap_semi_restriction semi direction wrong semi-direction. ");
    }
}

//
// interpolation
//

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_x (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic+1,jc,kc))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc+1,kc))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_z (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w2 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc,kc+1))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xy (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_yz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sig,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sig,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
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
            Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
            Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
            Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
            Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
            Real w5 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
            Real w6 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
            fine(i,j,k) += (w1 * aa_interp_face_yz(crse,sig,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * aa_interp_face_yz(crse,sig,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * aa_interp_face_xz(crse,sig,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * aa_interp_face_xz(crse,sig,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * aa_interp_face_xy(crse,sig,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * aa_interp_face_xy(crse,sig,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd and k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd and k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk, int idir) noexcept
{
    if (idir == 2 ) 
    {
    	if (!msk(i,j,k)) {
            int ic = amrex::coarsen(i,2);
            int jc = amrex::coarsen(j,2);
            int kc = k;
            bool i_is_odd = (ic*2 != i);
            bool j_is_odd = (jc*2 != j);

            if (i_is_odd and j_is_odd) {
                // Node on a X-Y face
                fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node 
                fine(i,j,k) += crse(ic,jc,kc);
            }
	}
    } else if (idir ==1 ){
        if (!msk(i,j,k)) {
            int ic = amrex::coarsen(i,2);
            int jc = j;
            int kc = amrex::coarsen(k,2);
            bool i_is_odd = (ic*2 != i);
            bool k_is_odd = (kc*2 != k);

            if (i_is_odd and k_is_odd) {
                // Node on a X-Z face
	        fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,kc);
            }
        }
    } else if (idir == 0 ) {
        if (!msk(i,j,k)) {
            int ic = i;
            int jc = amrex::coarsen(j,2);
            int kc = amrex::coarsen(k,2);
            bool j_is_odd = (jc*2 != j);
            bool k_is_odd = (kc*2 != k);

            if (j_is_odd and k_is_odd) {
                // Node on a Y-Z face
                fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,kc);
            }
        }
    } else { 
        amrex::Abort("mlndlap_semi_interpolation semi direction wrong semi-direction. ");
    }
}
namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xy (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        Real w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        Real w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xz (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        Real w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        Real w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_yz (Array4<Real const> const& crse,
                            Array4<Real const> const& sigy, Array4<Real const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        Real w2 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        Real w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        Real w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sigy,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sigy,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sigx,
                           Array4<Real const> const& sigy, Array4<Real const> const& sigz,
                           Array4<int const> const& msk) noexcept
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
            Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
            Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
            Real w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
            Real w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
            Real w5 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
            Real w6 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
            fine(i,j,k) += (w1 * ha_interp_face_yz(crse,sigy,sigz,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * ha_interp_face_yz(crse,sigy,sigz,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * ha_interp_face_xz(crse,sigx,sigz,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * ha_interp_face_xz(crse,sigx,sigz,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd and k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += ha_interp_face_yz(crse,sigy,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd and k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += ha_interp_face_xz(crse,sigx,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += ha_interp_face_xy(crse,sigx,sigy,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sigx,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sigy,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sigz,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

//
// rhs & u
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Box const& nodal_domain,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    Real facx = 0.25*dxinv[0];
    Real facy = 0.25*dxinv[1];
    Real facz = 0.25*dxinv[2];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {

        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;
        Real zero_klo = 1.0;
        Real zero_khi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0;
        if ((bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)
            and k == domlo.z) zero_klo = 0.0;
        if ((bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)
            and k == domhi.z) zero_khi = 0.0;

        rhs(i,j,k) = facx*(-vel(i-1,j-1,k-1,0)*zero_jlo*zero_klo+vel(i,j-1,k-1,0)*zero_jlo*zero_klo
                           -vel(i-1,j  ,k-1,0)*zero_jhi*zero_klo+vel(i,j  ,k-1,0)*zero_jhi*zero_klo
                           -vel(i-1,j-1,k  ,0)*zero_jlo*zero_khi+vel(i,j-1,k  ,0)*zero_jlo*zero_khi
                           -vel(i-1,j  ,k  ,0)*zero_jhi*zero_khi+vel(i,j  ,k  ,0)*zero_jhi*zero_khi)

                   + facy*(-vel(i-1,j-1,k-1,1)*zero_ilo*zero_klo-vel(i,j-1,k-1,1)*zero_ihi*zero_klo
                           +vel(i-1,j  ,k-1,1)*zero_ilo*zero_klo+vel(i,j  ,k-1,1)*zero_ihi*zero_klo
                           -vel(i-1,j-1,k  ,1)*zero_ilo*zero_khi-vel(i,j-1,k  ,1)*zero_ihi*zero_khi
                           +vel(i-1,j  ,k  ,1)*zero_ilo*zero_khi+vel(i,j  ,k  ,1)*zero_ihi*zero_khi)

                   + facz*(-vel(i-1,j-1,k-1,2)*zero_ilo*zero_jlo-vel(i,j-1,k-1,2)*zero_ihi*zero_jlo
                           -vel(i-1,j  ,k-1,2)*zero_ilo*zero_jhi-vel(i,j  ,k-1,2)*zero_ihi*zero_jhi
                           +vel(i-1,j-1,k  ,2)*zero_ilo*zero_jlo+vel(i,j-1,k  ,2)*zero_ihi*zero_jlo
                           +vel(i-1,j  ,k  ,2)*zero_ilo*zero_jhi+vel(i,j  ,k  ,2)*zero_ihi*zero_jhi);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc (int i, int j, int k, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    Real r;
    if (msk(i,j,k)) {
        r = 0.0;
    } else {
        r = 0.125 * (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) +
                     rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ));
    }
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.25*dxinv[0];
    Real facy = 0.25*dxinv[1];
    Real facz = 0.25*dxinv[2];
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_compute_fine_contrib (int i, int j, int k, Box const& fvbx,
                                        Array4<Real> const& frh, Array4<Real const> const& vel,
                                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                        Box const& nodal_domain,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    IntVect iv(i,j,k);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv))
    {
        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;
        Real zero_klo = 1.0;
        Real zero_khi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0;
        if ((bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)
            and k == domlo.z) zero_klo = 0.0;
        if ((bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)
            and k == domhi.z) zero_khi = 0.0;

        frh(i,j,k) = 0.25*dxinv[0]*(-vel(i-1,j-1,k-1,0)*zero_jlo*zero_klo+vel(i,j-1,k-1,0)*zero_jlo*zero_klo
                                    -vel(i-1,j  ,k-1,0)*zero_jhi*zero_klo+vel(i,j  ,k-1,0)*zero_jhi*zero_klo
                                    -vel(i-1,j-1,k  ,0)*zero_jlo*zero_khi+vel(i,j-1,k  ,0)*zero_jlo*zero_khi
                                    -vel(i-1,j  ,k  ,0)*zero_jhi*zero_khi+vel(i,j  ,k  ,0)*zero_jhi*zero_khi)
                   + 0.25*dxinv[1]*(-vel(i-1,j-1,k-1,1)*zero_ilo*zero_klo-vel(i,j-1,k-1,1)*zero_ihi*zero_klo
                                    +vel(i-1,j  ,k-1,1)*zero_ilo*zero_klo+vel(i,j  ,k-1,1)*zero_ihi*zero_klo
                                    -vel(i-1,j-1,k  ,1)*zero_ilo*zero_khi-vel(i,j-1,k  ,1)*zero_ihi*zero_khi
                                    +vel(i-1,j  ,k  ,1)*zero_ilo*zero_khi+vel(i,j  ,k  ,1)*zero_ihi*zero_khi)
                   + 0.25*dxinv[2]*(-vel(i-1,j-1,k-1,2)*zero_ilo*zero_jlo-vel(i,j-1,k-1,2)*zero_ihi*zero_jlo
                                    -vel(i-1,j  ,k-1,2)*zero_ilo*zero_jhi-vel(i,j  ,k-1,2)*zero_ihi*zero_jhi
                                    +vel(i-1,j-1,k  ,2)*zero_ilo*zero_jlo+vel(i,j-1,k  ,2)*zero_ihi*zero_jlo
                                    +vel(i-1,j  ,k  ,2)*zero_ilo*zero_jhi+vel(i,j  ,k  ,2)*zero_ihi*zero_jhi);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_add_fine_contrib (int i, int j, int k, Box const& fvbx,
                                    Array4<Real> const& rhs, Array4<Real const> const& frh,
                                    Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = Real(0.125);
    constexpr Real chip = Real(0.5);
    constexpr Real chip2 = Real(0.25);
    constexpr Real chip3 = Real(0.125);

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    IntVect iv(ii,jj,kk);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,kk))
    {
        rhs(i,j,k) +=
            rfd*(frh(ii,jj,kk)
                 + chip*(frh(ii,jj,kk-1)+frh(ii,jj,kk+1)
                         +frh(ii,jj-1,kk)+frh(ii,jj+1,kk)
                         +frh(ii-1,jj,kk)+frh(ii+1,jj,kk))
                 + chip2*(frh(ii,jj-1,kk-1)+frh(ii,jj+1,kk-1)+frh(ii,jj-1,kk+1)+frh(ii,jj+1,kk+1)
                          +frh(ii-1,jj,kk-1)+frh(ii+1,jj,kk-1)+frh(ii-1,jj,kk+1)+frh(ii+1,jj,kk+1)
                          +frh(ii-1,jj-1,kk)+frh(ii+1,jj-1,kk)+frh(ii-1,jj+1,kk)+frh(ii+1,jj+1,kk))
                 + chip3*(frh(ii-1,jj-1,kk-1)+frh(ii+1,jj-1,kk-1)
                          +frh(ii-1,jj+1,kk-1)+frh(ii+1,jj+1,kk-1)
                          +frh(ii-1,jj-1,kk+1)+frh(ii+1,jj-1,kk+1)
                          +frh(ii-1,jj+1,kk+1)+frh(ii+1,jj+1,kk+1)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int i, int j, int k, Box const& fvbx,
                                Array4<Real> const& rhs, Array4<Real const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    constexpr Real fac[] = {Real(0.125), Real(0.375), Real(0.375), Real(0.125)};
    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    IntVect iv(ii,jj,kk);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,kk))
    {
        Real r = 0.0;
        for (int koff = -2; koff <= 1; ++koff) {
        for (int joff = -2; joff <= 1; ++joff) {
        for (int ioff = -2; ioff <= 1; ++ioff) {
            r += cc(ii+ioff,jj+joff,kk+koff) * fac[ioff+2]*fac[joff+2]*fac[koff+2];
        }}}
        rhs(i,j,k) += r;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int i, int j, int k, Array4<Real> const& rhs,
                              Array4<Real const> const& vel, Array4<Real const> const& fc,
                              Array4<Real const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                              bool neumann_doubling) noexcept
{
    if (!dmsk(i,j,k) and ndmsk(i,j,k) == crse_fine_node) {
        Real facx = Real(0.25) * dxinv[0];
        Real facy = Real(0.25) * dxinv[1];
        Real facz = Real(0.25) * dxinv[2];
        Real r = fc(i,j,k);
        if (rhcc) {
            r += Real(0.125)*((Real(1.)-ccmsk(i-1,j-1,k-1)) * rhcc(i-1,j-1,k-1)
                            + (Real(1.)-ccmsk(i  ,j-1,k-1)) * rhcc(i  ,j-1,k-1)
                            + (Real(1.)-ccmsk(i-1,j  ,k-1)) * rhcc(i-1,j  ,k-1)
                            + (Real(1.)-ccmsk(i  ,j  ,k-1)) * rhcc(i  ,j  ,k-1)
                            + (Real(1.)-ccmsk(i-1,j-1,k  )) * rhcc(i-1,j-1,k  )
                            + (Real(1.)-ccmsk(i  ,j-1,k  )) * rhcc(i  ,j-1,k  )
                            + (Real(1.)-ccmsk(i-1,j  ,k  )) * rhcc(i-1,j  ,k  )
                            + (Real(1.)-ccmsk(i  ,j  ,k  )) * rhcc(i  ,j  ,k  ));
        }
        if (ccmsk(i-1,j-1,k-1) == crse_cell) {
            r += - facx*vel(i-1,j-1,k-1,0)
                 - facy*vel(i-1,j-1,k-1,1)
                 - facz*vel(i-1,j-1,k-1,2);
        }
        if (ccmsk(i,j-1,k-1) == crse_cell) {
            r += + facx*vel(i  ,j-1,k-1,0)
                 - facy*vel(i  ,j-1,k-1,1)
                 - facz*vel(i  ,j-1,k-1,2);
        }
        if (ccmsk(i-1,j,k-1) == crse_cell) {
            r += - facx*vel(i-1,j  ,k-1,0)
                 + facy*vel(i-1,j  ,k-1,1)
                 - facz*vel(i-1,j  ,k-1,2);
        }
        if (ccmsk(i,j,k-1) == crse_cell) {
            r += + facx*vel(i  ,j  ,k-1,0)
                 + facy*vel(i  ,j  ,k-1,1)
                 - facz*vel(i  ,j  ,k-1,2);
        }
        if (ccmsk(i-1,j-1,k) == crse_cell) {
            r += - facx*vel(i-1,j-1,k  ,0)
                 - facy*vel(i-1,j-1,k  ,1)
                 + facz*vel(i-1,j-1,k  ,2);
        }
        if (ccmsk(i,j-1,k) == crse_cell) {
            r += + facx*vel(i  ,j-1,k  ,0)
                 - facy*vel(i  ,j-1,k  ,1)
                 + facz*vel(i  ,j-1,k  ,2);
        }
        if (ccmsk(i-1,j,k) == crse_cell) {
            r += - facx*vel(i-1,j  ,k  ,0)
                 + facy*vel(i-1,j  ,k  ,1)
                 + facz*vel(i-1,j  ,k  ,2);
        }
        if (ccmsk(i,j,k) == crse_cell) {
            r += + facx*vel(i  ,j  ,k  ,0)
                 + facy*vel(i  ,j  ,k  ,1)
                 + facz*vel(i  ,j  ,k  ,2);
        }

        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                                  bclo[0] == LinOpBCType::inflow)) {
                r *= Real(2.);
            } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                        bchi[0] == LinOpBCType::inflow)) {
                r *= Real(2.);
            }

            if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                                  bclo[1] == LinOpBCType::inflow)) {
                r *= Real(2.);
            } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                         bchi[1] == LinOpBCType::inflow)) {
                r *= Real(2.);
            }

            if (k == ndlo.z and ( bclo[2] == LinOpBCType::Neumann or
                                  bclo[2] == LinOpBCType::inflow)) {
                r *= Real(2.);
            } else if (k == ndhi.z and ( bchi[2] == LinOpBCType::Neumann or
                                         bchi[2] == LinOpBCType::inflow)) {
                r *= Real(2.);
            }
        }

        rhs(i,j,k) = r;
    }
}

//
// residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int i, int j, int k, Array4<Real> const& resid,
                         Array4<Real const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                         bool neumann_doubling) noexcept
{
    if ((msk(i-1,j-1,k-1) == 0 or
         msk(i  ,j-1,k-1) == 0 or
         msk(i-1,j  ,k-1) == 0 or
         msk(i  ,j  ,k-1) == 0 or
         msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0) and
        (msk(i-1,j-1,k-1) == 0 or
         msk(i  ,j-1,k-1) == 0 or
         msk(i-1,j  ,k-1) == 0 or
         msk(i  ,j  ,k-1) == 0 or
         msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0))
    {
        Real fac = Real(1.0);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                                  bclo[0] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                        bchi[0] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            }
            if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                                  bclo[1] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                         bchi[1] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            }
            if (k == ndlo.z and ( bclo[2] == LinOpBCType::Neumann or
                                  bclo[2] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            } else if (k == ndhi.z and ( bchi[2] == LinOpBCType::Neumann or
                                         bchi[2] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            }
        }
        resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac;
    } else {
        resid(i,j,k) = Real(0.);
    }
}

//
// sync residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_Ax (int i, int j, int k, Box const& fvbx, Array4<Real> const& Ax,
                          Array4<Real const> const& x, Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j,k);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv)) {
        Real facx = Real(1./36.)*dxinv[0]*dxinv[0];
        Real facy = Real(1./36.)*dxinv[1]*dxinv[1];
        Real facz = Real(1./36.)*dxinv[2]*dxinv[2];
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + Real(2.)*facy + Real(2.)*facz;
        Real f2xmy2z = Real(2.)*facx - facy + Real(2.)*facz;
        Real f2x2ymz = Real(2.)*facx + Real(2.)*facy - facz;
        Real f4xm2ym2z = Real(4.)*facx - Real(2.)*facy - Real(2.)*facz;
        Real fm2x4ym2z = -Real(2.)*facx + Real(4.)*facy - Real(2.)*facz;
        Real fm2xm2y4z = -Real(2.)*facx - Real(2.)*facy + Real(4.)*facz;
        Ax(i,j,k) = x(i,j,k)*Real(-4.)*fxyz*
            (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
            //
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
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_contrib (int i, int j, int k, Array4<Real> const& f,
                               Array4<Real const> const& Ax, Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = Real(0.125);
    constexpr Real chip = Real(0.5);
    constexpr Real chip2 = Real(0.25);
    constexpr Real chip3 = Real(0.125);

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    if (msk(ii,jj,kk)) {
        f(i,j,k) += rfd*
            (Ax(ii,jj,kk)
             + chip*(Ax(ii,jj,kk-1)+Ax(ii,jj,kk+1)
                     +Ax(ii,jj-1,kk)+Ax(ii,jj+1,kk)
                     +Ax(ii-1,jj,kk)+Ax(ii+1,jj,kk))
             + chip2*(Ax(ii,jj-1,kk-1)+Ax(ii,jj+1,kk-1)+Ax(ii,jj-1,kk+1)+Ax(ii,jj+1,kk+1)
                      +Ax(ii-1,jj,kk-1)+Ax(ii+1,jj,kk-1)+Ax(ii-1,jj,kk+1)+Ax(ii+1,jj,kk+1)
                      +Ax(ii-1,jj-1,kk)+Ax(ii+1,jj-1,kk)+Ax(ii-1,jj+1,kk)+Ax(ii+1,jj+1,kk))
             + chip3*(Ax(ii-1,jj-1,kk-1)+Ax(ii+1,jj-1,kk-1)
                      +Ax(ii-1,jj+1,kk-1)+Ax(ii+1,jj+1,kk-1)
                      +Ax(ii-1,jj-1,kk+1)+Ax(ii+1,jj-1,kk+1)
                      +Ax(ii-1,jj+1,kk+1)+Ax(ii+1,jj+1,kk+1)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int i, int j, int k, Array4<Real> const& res,
                             Array4<Real const> const& phi, Array4<Real const> const& rhs,
                             Array4<Real const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<Real const> const& fc,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                             bool neumann_doubling) noexcept
{
    if (!dmsk(i,j,k) and ndmsk(i,j,k) == crse_fine_node) {

        Real facx = Real(1./36.)*dxinv[0]*dxinv[0];
        Real facy = Real(1./36.)*dxinv[1]*dxinv[1];
        Real facz = Real(1./36.)*dxinv[2]*dxinv[2];

        Real Ax = Real(0.);
        if (ccmsk(i-1,j-1,k-1) == crse_cell) {
            Ax += sig(i-1,j-1,k-1)*(facx*(Real(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                         +Real(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                         +Real(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                         +         (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                  + facy*(Real(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                         +Real(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                         +Real(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                         +         (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1)))
                                  + facz*(Real(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                         +Real(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                         +Real(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                         +         (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  ))));
        }
        if (ccmsk(i,j-1,k-1) == crse_cell) {
            Ax += sig(i,j-1,k-1)*(facx*(Real(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +Real(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                + facy*(Real(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                       +Real(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1)))
                                + facz*(Real(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                       +Real(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                       +         (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  ))));
        }
        if (ccmsk(i-1,j,k-1) == crse_cell) {
            Ax += sig(i-1,j,k-1)*(facx*(Real(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                       +Real(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1)))
                                 + facy*(Real(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                        +Real(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                        +Real(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                        +         (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1)))
                                 + facz*(Real(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                        +Real(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                        +Real(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                        +         (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  ))));
        }
        if (ccmsk(i,j,k-1) == crse_cell) {
            Ax += sig(i,j,k-1)*(facx*(Real(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +Real(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1)))
                              + facy*(Real(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1)))
                              + facz*(Real(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                     +         (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  ))));
        }
        if (ccmsk(i-1,j-1,k) == crse_cell) {
            Ax += sig(i-1,j-1,k)*(facx*(Real(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +Real(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                       +         (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1)))
                                + facy*(Real(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                       +Real(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                       +         (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1)))
                                + facz*(Real(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                       +Real(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                       +Real(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                       +         (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  ))));
        }
        if (ccmsk(i,j-1,k) == crse_cell) {
            Ax += sig(i,j-1,k)*(facx*(Real(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                     +Real(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1)))
                              + facy*(Real(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1)))
                              + facz*(Real(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                     +         (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  ))));
        }
        if (ccmsk(i-1,j,k) == crse_cell) {
            Ax += sig(i-1,j,k)*(facx*(Real(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +Real(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1)))
                              + facy*(Real(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1)))
                              + facz*(Real(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +Real(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                     +Real(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                     +         (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  ))));
        }
        if (ccmsk(i,j,k) == crse_cell) {
            Ax += sig(i,j,k)*(facx*(Real(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +Real(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                   +Real(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1)))
                            + facy*(Real(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                   +Real(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                   +Real(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1)))
                            + facz*(Real(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                   +Real(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                   +Real(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                   +         (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  ))));
        }

        Real Axf = fc(i,j,k);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);

            if (i == ndlo.x and (bclo[0] == LinOpBCType::Neumann or
                                 bclo[0] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            } else if (i== ndhi.x and (bchi[0] == LinOpBCType::Neumann or
                                       bchi[0] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            }

            if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                                  bclo[1] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            } else if (j == ndhi.y and (bchi[1] == LinOpBCType::Neumann or
                                        bchi[1] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            }

            if (k == ndlo.z and (bclo[2] == LinOpBCType::Neumann or
                                 bclo[2] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            } else if (k == ndhi.z and (bchi[2] == LinOpBCType::Neumann or
                                        bchi[2] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            }
        }

        res(i,j,k) = rhs(i,j,k) - (Ax + Axf);
    }
}

//
// RAP
//

namespace {

    constexpr int ist_000 = 0;
    constexpr int ist_p00 = 1;
    constexpr int ist_0p0 = 2;
    constexpr int ist_00p = 3;
    constexpr int ist_pp0 = 4;
    constexpr int ist_p0p = 5;
    constexpr int ist_0pp = 6;
    constexpr int ist_ppp = 7;
    constexpr int ist_inv = 8;
    constexpr int n_sten  = 9;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const& bx, Array4<Real> const& sten,
                          Array4<Real const> const& sig,
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
        // i+1,j,k
        sten(i,j,k,ist_p00) = f4xm2ym2z * (sig(i,j-1,k-1)+sig(i,j,k-1)+sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k: sten(i-1,j,k,ist_p00)

        // i,j+1,k
        sten(i,j,k,ist_0p0) = fm2x4ym2z * (sig(i-1,j,k-1)+sig(i,j,k-1)+sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k: sten(i,j-1,k,ist_0p0)

        // i,j,k+1
        sten(i,j,k,ist_00p) = fm2xm2y4z * (sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        // i,j,k-1: sten(i,j,k-1,ist_00p)

        // i+1,j+1,k
        sten(i,j,k,ist_pp0) = f2x2ymz * (sig(i,j,k-1)+sig(i,j,k));
        // i-1,j-1,k: sten(i-1,j-1,k,ist_pp0)
        // i+1,j-1,k: sten(i  ,j-1,k,ist_pp0)
        // i-1,j+1,k: sten(i-1,j  ,k,ist_pp0)

        // i+1,j,k+1
        sten(i,j,k,ist_p0p) = f2xmy2z * (sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k-1: sten(i-1,j,k-1,ist_p0p)
        // i+1,j,k-1: sten(i  ,j,k-1,ist_p0p)
        // i-1,j,k+1: sten(i-1,j,k  ,ist_p0p)

        // i,j+1,k+1
        sten(i,j,k,ist_0pp) = fmx2y2z * (sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k-1: sten(i,j-1,k-1,ist_0pp)
        // i,j+1,k-1: sten(i,j  ,k-1,ist_0pp)
        // i,j-1,k+1: sten(i,j-1,k  ,ist_0pp)

        // i+1,j+1,k+1
        sten(i,j,k,ist_ppp) = fxyz * sig(i,j,k);
        // i-1,j-1,k-1: sten(i-1,j-1,k-1,ist_ppp)
        // i+1,j-1,k-1: sten(i  ,j-1,k-1,ist_ppp)
        // i-1,j+1,k-1: sten(i-1,j  ,k-1,ist_ppp)
        // i+1,j+1,k-1: sten(i  ,j  ,k-1,ist_ppp)
        // i-1,j-1,k+1: sten(i-1,j-1,k  ,ist_ppp)
        // i+1,j-1,k+1: sten(i  ,j-1,k  ,ist_ppp)
        // i-1,j+1,k+1: sten(i-1,j  ,k  ,ist_ppp)
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int i, int j, int k, Array4<Real> const& sten) noexcept
{
    sten(i,j,k,ist_000) = -(sten(i-1,j,k,ist_p00) + sten(i,j,k,ist_p00)
                          + sten(i,j-1,k,ist_0p0) + sten(i,j,k,ist_0p0)
                          + sten(i,j,k-1,ist_00p) + sten(i,j,k,ist_00p)
                          + sten(i-1,j-1,k,ist_pp0) + sten(i,j-1,k,ist_pp0)
                          + sten(i-1,j,k,ist_pp0) + sten(i,j,k,ist_pp0)
                          + sten(i-1,j,k-1,ist_p0p) + sten(i,j,k-1,ist_p0p)
                          + sten(i-1,j,k,ist_p0p) + sten(i,j,k,ist_p0p)
                          + sten(i,j-1,k-1,ist_0pp) + sten(i,j,k-1,ist_0pp)
                          + sten(i,j-1,k,ist_0pp) + sten(i,j,k,ist_0pp)
                          + sten(i-1,j-1,k-1,ist_ppp) + sten(i,j-1,k-1,ist_ppp)
                          + sten(i-1,j,k-1,ist_ppp) + sten(i,j,k-1,ist_ppp)
                          + sten(i-1,j-1,k,ist_ppp) + sten(i,j-1,k,ist_ppp)
                          + sten(i-1,j,k,ist_ppp) + sten(i,j,k,ist_ppp));
    sten(i,j,k,ist_inv) = 1.0 /
        (  amrex::Math::abs(sten(i-1,j,k,ist_p00)) + amrex::Math::abs(sten(i,j,k,ist_p00))
         + amrex::Math::abs(sten(i,j-1,k,ist_0p0)) + amrex::Math::abs(sten(i,j,k,ist_0p0))
         + amrex::Math::abs(sten(i,j,k-1,ist_00p)) + amrex::Math::abs(sten(i,j,k,ist_00p))
         + amrex::Math::abs(sten(i-1,j-1,k,ist_pp0)) + amrex::Math::abs(sten(i,j-1,k,ist_pp0))
         + amrex::Math::abs(sten(i-1,j,k,ist_pp0)) + amrex::Math::abs(sten(i,j,k,ist_pp0))
         + amrex::Math::abs(sten(i-1,j,k-1,ist_p0p)) + amrex::Math::abs(sten(i,j,k-1,ist_p0p))
         + amrex::Math::abs(sten(i-1,j,k,ist_p0p)) + amrex::Math::abs(sten(i,j,k,ist_p0p))
         + amrex::Math::abs(sten(i,j-1,k-1,ist_0pp)) + amrex::Math::abs(sten(i,j,k-1,ist_0pp))
         + amrex::Math::abs(sten(i,j-1,k,ist_0pp)) + amrex::Math::abs(sten(i,j,k,ist_0pp))
         + amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp)) + amrex::Math::abs(sten(i,j-1,k-1,ist_ppp))
         + amrex::Math::abs(sten(i-1,j,k-1,ist_ppp)) + amrex::Math::abs(sten(i,j,k-1,ist_ppp))
         + amrex::Math::abs(sten(i-1,j-1,k,ist_ppp)) + amrex::Math::abs(sten(i,j-1,k,ist_ppp))
         + amrex::Math::abs(sten(i-1,j,k,ist_ppp)) + amrex::Math::abs(sten(i,j,k,ist_ppp)) + eps);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int i, int j, int k, Array4<Real> const& csten,
                          Array4<Real const> const& fsten) noexcept
{
    auto interp_from_mmm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };
    amrex::ignore_unused(interp_from_mmm_to);

    auto interp_from_pmm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mpm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mmp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_pmp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_p00)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) * fsten(i_,j_  ,k_,ist_inv);
        return p;
    };

    auto interp_from_mpp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = Real(1.);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p *= amrex::Math::abs(fsten(i_  ,j_  ,k_,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_0mm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0mp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0pm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0pp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_m0m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_p0m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_m0p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_p0p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_mm0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
      Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                            +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_mp0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_pm0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_pp0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_00m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_00p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_0m0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_0p0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_m00_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_p00_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = amrex::Math::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        Real w2 = amrex::Math::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto Ammm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_-1,ist_ppp);
    };
    amrex::ignore_unused(Ammm);

    auto A0mm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_-1,ist_0pp);
    };

    auto Apmm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_-1,ist_ppp);
    };

    auto Am0m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_-1,ist_p0p);
    };

    auto A00m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_00p);
    };

    auto Ap0m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_p0p);
    };

    auto Ampm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_-1,ist_ppp);
    };

    auto A0pm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_0pp);
    };

    auto Appm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_ppp);
    };

    auto Amm0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_  ,ist_pp0);
    };

    auto A0m0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_0p0);
    };

    auto Apm0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_pp0);
    };

    auto Am00 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_p00);
    };

    auto A000 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_000);
    };

    auto Ap00 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_p00);
    };

    auto Amp0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_pp0);
    };

    auto A0p0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_0p0);
    };

    auto App0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_pp0);
    };

    auto Ammp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_  ,ist_ppp);
    };

    auto A0mp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_0pp);
    };

    auto Apmp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_ppp);
    };

    auto Am0p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_p0p);
    };

    auto A00p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_00p);
    };

    auto Ap0p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_p0p);
    };

    auto Ampp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_ppp);
    };

    auto A0pp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_0pp);
    };

    auto Appp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_ppp);
    };

    auto restrict_from_mmm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-2,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) * fsten(ii_-1,jj_-1,kk_-1,ist_inv);
        return r;
    };
    amrex::ignore_unused(restrict_from_mmm_to);

    auto restrict_from_0mm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_,jj_-2,kk_-1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_,jj_-1,kk_-1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_,jj_-1,kk_-2,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_,jj_-1,kk_-1,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_pmm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) * fsten(ii_+1,jj_-1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_m0m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_-2,jj_,kk_-1,ist_p00)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_-1,jj_,kk_-1,ist_p00)) / (amrex::Math::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_-1,jj_,kk_-2,ist_00p)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_-1,jj_,kk_-1,ist_00p)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_00m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_,jj_,kk_-2,ist_00p));
        Real w2 = amrex::Math::abs(fsten(ii_,jj_,kk_-1,ist_00p));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_p0m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_  ,jj_,kk_-1,ist_p00)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_+1,jj_,kk_-1,ist_p00)) / (amrex::Math::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_+1,jj_,kk_-2,ist_00p)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_+1,jj_,kk_-1,ist_00p)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mpm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_+1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) * fsten(ii_-1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_0pm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_,jj_  ,kk_-1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_,jj_+1,kk_-1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_,jj_+1,kk_-2,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_,jj_+1,kk_-1,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_ppm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_+1,kk_-1,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) * fsten(ii_+1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_mm0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_-2,jj_-1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_-1,jj_-1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_-1,jj_-2,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_-1,jj_-1,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_0m0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_,jj_-2,kk_,ist_0p0));
        Real w2 = amrex::Math::abs(fsten(ii_,jj_-1,kk_,ist_0p0));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_pm0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_  ,jj_-1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_+1,jj_-1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_+1,jj_-2,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_+1,jj_-1,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_m00_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_-2,jj_,kk_,ist_p00));
        Real w2 = amrex::Math::abs(fsten(ii_-1,jj_,kk_,ist_p00));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_000_to = [] (int /*ii_*/, int /*jj_*/, int /*kk_*/) -> Real {
        return Real(1.);
    };

    auto restrict_from_p00_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_  ,jj_,kk_,ist_p00));
        Real w2 = amrex::Math::abs(fsten(ii_+1,jj_,kk_,ist_p00));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_mp0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_-2,jj_+1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_-1,jj_+1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_-1,jj_  ,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_-1,jj_+1,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_0p0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_,jj_  ,kk_,ist_0p0));
        Real w2 = amrex::Math::abs(fsten(ii_,jj_+1,kk_,ist_0p0));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_pp0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_  ,jj_+1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_+1,jj_+1,kk_,ist_p00)) / (amrex::Math::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_+1,jj_  ,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_+1,jj_+1,kk_,ist_0p0)) / (amrex::Math::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mmp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-2,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_-1,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) * fsten(ii_-1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0mp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_,jj_-2,kk_+1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_,jj_-1,kk_+1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_,jj_-1,kk_  ,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_,jj_-1,kk_+1,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_pmp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_+1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-2,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_+1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_+1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_-1,kk_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) * fsten(ii_+1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_m0p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_-2,jj_,kk_+1,ist_p00)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_-1,jj_,kk_+1,ist_p00)) / (amrex::Math::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_-1,jj_,kk_  ,ist_00p)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_-1,jj_,kk_+1,ist_00p)) / (amrex::Math::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_00p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = amrex::Math::abs(fsten(ii_,jj_,kk_  ,ist_00p));
        Real w2 = amrex::Math::abs(fsten(ii_,jj_,kk_+1,ist_00p));
        if (w1 == Real(0.) and w2 == Real(0.)) {
            return Real(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_p0p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_  ,jj_,kk_+1,ist_p00)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_+1,jj_,kk_+1,ist_p00)) / (amrex::Math::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_+1,jj_,kk_  ,ist_00p)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_+1,jj_,kk_+1,ist_00p)) / (amrex::Math::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p))
                                                              +amrex::Math::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mpp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_+1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_+1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_  ,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-2,jj_+1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_+1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_+1,kk_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) * fsten(ii_-1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0pp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = amrex::Math::abs(fsten(ii_,jj_  ,kk_+1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) + eps);
        Real w1p = amrex::Math::abs(fsten(ii_,jj_+1,kk_+1,ist_0p0)) / (amrex::Math::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps);
        Real w2m = amrex::Math::abs(fsten(ii_,jj_+1,kk_  ,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) + eps);
        Real w2p = amrex::Math::abs(fsten(ii_,jj_+1,kk_+1,ist_00p)) / (amrex::Math::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp))
                                                              +amrex::Math::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps);
        Real wmm = amrex::Math::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp)) * (Real(1.) + w1m + w2m);
        Real wpm = amrex::Math::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) * (Real(1.) + w1p + w2m);
        Real wmp = amrex::Math::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) * (Real(1.) + w1m + w2p);
        Real wpp = amrex::Math::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) * (Real(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_ppp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = Real(1.);
        r += amrex::Math::abs(fsten(ii_  ,jj_+1,kk_+1,ist_p00)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_  ,kk_+1,ist_0p0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_+1,kk_  ,ist_00p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_  ,kk_+1,ist_pp0)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_  ,jj_+1,kk_  ,ist_p0p)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp)) + eps);
        r += amrex::Math::abs(fsten(ii_+1,jj_  ,kk_  ,ist_0pp)) /
           ( amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + amrex::Math::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp)) + eps);
        r *= amrex::Math::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp)) * fsten(ii_+1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    Array3D<Real,-1,1,-1,1,-1,1> p;
    Array3D<Real,-1,1,-1,1,-1,1> ap;
    Real cs1, cs2, cs3, cs4;

    // csten(i,j,k,ist_p00)
    int iii = ii;
    int jjj = jj;
    int kkk = kk;;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk-1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk-1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk  );
    p( 0, 0, 0) = Real(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk  );
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,+1) = interp_from_pmm_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,+1) = interp_from_0mm_to(iii+2,jjj+1,kkk+1);
    ap(0,-1,-1) =
                     Ap00(iii,jjj-1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk-1) * p(-1, 0, 0);
    ap(1,-1,-1) =
                     A000(iii+1,jjj-1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk-1) * p( 0, 0, 0);
    ap(0,0,-1) =
                     Apm0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk-1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk-1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk-1) * p(-1,+1, 0);
    ap(1,0,-1) =
                     A0m0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,+1, 0);
    ap(0,1,-1) =
                     Apm0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,+1, 0);
    ap(1,1,-1) =
                     A0m0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,+1, 0);
    ap(0,-1,0) =
                     Ap0m(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj-1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,+1);
    ap(1,-1,0) =
                     A00m(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj-1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj-1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj-1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,+1);
    ap(0,0,0) =
                     Apmm(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0m(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appm(iii,jjj,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Ap00(iii,jjj,kkk) * p(-1, 0, 0)
      +              App0(iii,jjj,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj,kkk) * p(-1,-1,+1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,+1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,+1);
    ap(1,0,0) =
                     A0mm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0m(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              Apm0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj,kkk) * p(-1, 0, 0)
      +              Ap00(iii+1,jjj,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj,kkk) * p(-1,+1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,+1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,+1);
    ap(0,1,0) =
                     Apmm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,+1);
    ap(1,1,0) =
                     A0mm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              Apm0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk) * p(-1,+1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,+1);
    ap(0,-1,1) =
                     Ap0m(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,+1);
    ap(1,-1,1) =
                     A00m(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,+1);
    ap(0,0,1) =
                     Apmm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appm(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,+1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,+1);
    ap(1,0,1) =
                     A0mm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,+1);
    ap(0,1,1) =
                     Apmm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,+1);
    ap(1,1,1) =
                     A0mm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,+1);
    csten(i,j,k,ist_p00) = Real(0.125) *
      ( restrict_from_0mm_to(iii,jjj,kkk) * ap( 0,-1,-1)
      + restrict_from_pmm_to(iii,jjj,kkk) * ap(+1,-1,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_0p0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk  );
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii+1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     A0p0(iii-1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk-1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk-1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(1,0,-1) =
                     Amp0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     A000(iii-1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk-1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk-1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk-1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk-1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(1,1,-1) =
                     Am00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     A0pm(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pm(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(-1,-1, 0)
      +              A0p0(iii,jjj,kkk) * p( 0,-1, 0)
      +              App0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(-1,-1,+1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,+1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(1,0,0) =
                     Ampm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     A00m(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00m(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pm(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              A000(iii,jjj+1,kkk) * p( 0,-1, 0)
      +              Ap00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              A0p0(iii,jjj+1,kkk) * p( 0, 0, 0)
      +              App0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,+1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(1,1,0) =
                     Am0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     A0pm(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pm(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,+1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(1,0,1) =
                     Ampm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     A00m(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00m(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pm(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,+1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(1,1,1) =
                     Am0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,+1);
    csten(i,j,k,ist_0p0) = Real(0.125) *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_00p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj  ,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj  ,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj+1,kkk+2);
    ap(-1,-1,0) =
                     A00p(iii-1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(1,-1,0) =
                     Am0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              A00p(iii,jjj,kkk) * p( 0, 0,-1)
      +              Ap0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     A0mp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(1,1,0) =
                     Ammp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A000(iii-1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(1,-1,1) =
                     Am00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              A000(iii,jjj,kkk+1) * p( 0, 0,-1)
      +              Ap00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              A00p(iii,jjj,kkk+1) * p( 0, 0, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     A0m0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(1,1,1) =
                     Amm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,+1, 0);
    csten(i,j,k,ist_00p) = Real(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_pp0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk  );
    p( 0, 0, 0) = Real(1.);
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj+2,kkk+1);
    ap(0,0,-1) =
                     App0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk-1) * p(-1,-1, 0);
    ap(1,0,-1) =
                     A0p0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,1,-1) =
                     Ap00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(-1, 0, 0);
    ap(1,1,-1) =
                     A000(iii+1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,0,0) =
                     Appm(iii,jjj,kkk) * p(-1,-1,-1)
      +              App0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Appp(iii,jjj,kkk) * p(-1,-1,+1);
    ap(1,0,0) =
                     A0pm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,+1);
    ap(0,1,0) =
                     Ap0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,+1);
    ap(1,1,0) =
                     A00m(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,0,1) =
                     Appm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              App0(iii,jjj,kkk+1) * p(-1,-1,+1);
    ap(1,0,1) =
                     A0pm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,1,1) =
                     Ap0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,+1);
    ap(1,1,1) =
                     A00m(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1);
    cs1 = Real(0.125) *
      ( restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mp0)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk-1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk  );
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk  );
    p( 0,-1,+1) = interp_from_0pm_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii-1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     Amp0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk-1) * p(+1,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     Am00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     Ampm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     Am0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii-1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     Ampm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     Am0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    cs2 = Real(0.125) *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_pp0) = Real(0.5)*(cs1 + cs2);

    // csten(i,j,k,ist_p0p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk+2);
    ap(0,-1,0) =
                     Ap0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,-1);
    ap(1,-1,0) =
                     A00p(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Apmp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,-1);
    ap(1,1,0) =
                     A0mp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Ap00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(-1, 0, 0);
    ap(1,-1,1) =
                     A000(iii+1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Apm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,+1, 0);
    ap(1,1,1) =
                     A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0);
    cs1 = Real(0.125) *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_m0p)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj  ,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj-1,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj  ,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj+1,kkk+2);

    ap(-1,-1,0) =
                     Am0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(+1, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     Ammp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Am00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     Amm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    cs2 = Real(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_p0p) = Real(0.5)*(cs1+cs2);

    // csten(i,j,k,ist_0pp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     A0pp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(1,0,0) =
                     Ampp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     A00p(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(1,1,0) =
                     Am0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     A0p0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(1,0,1) =
                     Amp0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     A000(iii-1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(1,1,1) =
                     Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0);
    cs1 = Real(0.125) *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i,j+1,k,ist_0mp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj-2,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj-2,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     A0mp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(1,-1,0) =
                     Ammp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A0m0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(1,-1,1) =
                     Amm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    cs2 = Real(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    csten(i,j,k,ist_0pp) = Real(0.5)*(cs1+cs2);

    // csten(i,j,k,ist_ppp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk+2);
    p( 0, 0, 0) = Real(1.);
    ap(0,0,0) =
                     Appp(iii,jjj,kkk) * p(-1,-1,-1);
    ap(1,0,0) =
                     A0pp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,-1);
    ap(0,1,0) =
                     Ap0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,-1);
    ap(1,1,0) =
                     A00p(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,0,1) =
                     App0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk+1) * p(-1,-1, 0);
    ap(1,0,1) =
                     A0p0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,1,1) =
                     Ap00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(-1, 0, 0);
    ap(1,1,1) =
                     A000(iii+1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0);
    cs1 = Real(0.125) *
      ( restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mpp)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     Ampp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     Am0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     Amp0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    cs2 = Real(0.125) *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    // alternative: csten(i,j+1,k,ist_pmp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj-2,kkk+2);
    p( 0, 0, 0) = Real(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj-1,kkk+2);
    ap(0,-1,0) =
                     Apmp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,+1,-1);
    ap(1,-1,0) =
                     A0mp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Apm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,+1, 0);
    ap(1,-1,1) =
                     A0m0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    cs3 = Real(0.125) *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    // alternative: csten(i+1,j+1,k,ist_mmp)
    iii = ii+2;
    jjj = jj+2;
    kkk = kk;
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj-2,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0, 0) = Real(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj-2,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     Ammp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Amm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,+1, 0);
    cs4 = Real(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1));

    csten(i,j,k,ist_ppp) = Real(0.25)*(cs1+cs2+cs3+cs4);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_sten (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        y(i,j,k) = x(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000)
            //
            +      x(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00)
            +      x(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00)
            //
            +      x(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0)
            +      x(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0)
            //
            +      x(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p)
            +      x(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p)
            //
            +      x(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0)
            +      x(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0)
            +      x(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0)
            +      x(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0)
            //
            +      x(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p)
            +      x(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p)
            +      x(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p)
            +      x(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p)
            //
            +      x(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp)
            +      x(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp)
            +      x(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp)
            +      x(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp)
            //
            +      x(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp)
            +      x(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp)
            +      x(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp)
            +      x(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp)
            +      x(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp)
            +      x(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp)
            +      x(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp)
            +      x(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (Box const& bx, Array4<Real> const& sol,
                                Array4<Real const> const& rhs,
                                Array4<Real const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else if (sten(i,j,k,ist_000) != 0.0) {
            Real Ax  = sol(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000)
                //
                +      sol(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00)
                +      sol(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00)
                //
                +      sol(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0)
                +      sol(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0)
                //
                +      sol(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p)
                +      sol(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p)
                //
                +      sol(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0)
                +      sol(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0)
                +      sol(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0)
                +      sol(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0)
                //
                +      sol(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p)
                +      sol(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p)
                +      sol(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p)
                +      sol(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p)
                //
                +      sol(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp)
                +      sol(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp)
                +      sol(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp)
                +      sol(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp)
                //
                +      sol(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp)
                +      sol(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp)
                +      sol(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp)
                +      sol(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp)
                +      sol(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp)
                +      sol(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp)
                +      sol(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp)
                +      sol(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp);

            sol(i,j,k) += (rhs(i,j,k) - Ax) / sten(i,j,k,ist_000);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int i, int j, int k, Array4<Real> const& fine,
                            Array4<Real const> const& crse, Array4<Real const> const& sten,
                            Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k) and sten(i,j,k,ist_000) != 0.0) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool ieven = ic*2 == i;
        bool jeven = jc*2 == j;
        bool keven = kc*2 == k;
        Real fv;
        if (ieven and jeven and keven) {
            fv = crse(ic,jc,kc);
        } else if (ieven and jeven) {
            Real w1 = amrex::Math::abs(sten(i,j,k-1,ist_00p));
            Real w2 = amrex::Math::abs(sten(i,j,k  ,ist_00p));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc,kc+1));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc,kc+1)) / (w1+w2);
            }
        } else if (ieven and keven) {
            Real w1 = amrex::Math::abs(sten(i,j-1,k,ist_0p0));
            Real w2 = amrex::Math::abs(sten(i,j  ,k,ist_0p0));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc+1,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc+1,kc)) / (w1+w2);
            }
        } else if (jeven and keven) {
            Real w1 = amrex::Math::abs(sten(i-1,j,k,ist_p00));
            Real w2 = amrex::Math::abs(sten(i  ,j,k,ist_p00));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic+1,jc,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic+1,jc,kc)) / (w1+w2);
            }
        } else if (ieven) {
            Real w1m = amrex::Math::abs(sten(i,j-1,k,ist_0p0)) / (amrex::Math::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +amrex::Math::abs(sten(i,j-1,k  ,ist_0pp)) + eps);
            Real w1p = amrex::Math::abs(sten(i,j  ,k,ist_0p0)) / (amrex::Math::abs(sten(i,j  ,k-1,ist_0pp))
                                                         +amrex::Math::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real w2m = amrex::Math::abs(sten(i,j,k-1,ist_00p)) / (amrex::Math::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +amrex::Math::abs(sten(i,j  ,k-1,ist_0pp)) + eps);
            Real w2p = amrex::Math::abs(sten(i,j,k  ,ist_00p)) / (amrex::Math::abs(sten(i,j-1,k  ,ist_0pp))
                                                         +amrex::Math::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real wmm = amrex::Math::abs(sten(i,j-1,k-1,ist_0pp)) * (1.0 + w1m + w2m);
            Real wpm = amrex::Math::abs(sten(i,j  ,k-1,ist_0pp)) * (1.0 + w1p + w2m);
            Real wmp = amrex::Math::abs(sten(i,j-1,k  ,ist_0pp)) * (1.0 + w1m + w2p);
            Real wpp = amrex::Math::abs(sten(i,j  ,k  ,ist_0pp)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic,jc+1,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic,jc+1,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (jeven) {
            Real w1m = amrex::Math::abs(sten(i-1,j,k,ist_p00)) / (amrex::Math::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +amrex::Math::abs(sten(i-1,j,k  ,ist_p0p)) + eps);
            Real w1p = amrex::Math::abs(sten(i  ,j,k,ist_p00)) / (amrex::Math::abs(sten(i  ,j,k-1,ist_p0p))
                                                         +amrex::Math::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real w2m = amrex::Math::abs(sten(i,j,k-1,ist_00p)) / (amrex::Math::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +amrex::Math::abs(sten(i  ,j,k-1,ist_p0p)) + eps);
            Real w2p = amrex::Math::abs(sten(i,j,k  ,ist_00p)) / (amrex::Math::abs(sten(i-1,j,k  ,ist_p0p))
                                                         +amrex::Math::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real wmm = amrex::Math::abs(sten(i-1,j,k-1,ist_p0p)) * (1.0 + w1m + w2m);
            Real wpm = amrex::Math::abs(sten(i  ,j,k-1,ist_p0p)) * (1.0 + w1p + w2m);
            Real wmp = amrex::Math::abs(sten(i-1,j,k  ,ist_p0p)) * (1.0 + w1m + w2p);
            Real wpp = amrex::Math::abs(sten(i  ,j,k  ,ist_p0p)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic+1,jc,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (keven) {
            Real w1m = amrex::Math::abs(sten(i-1,j,k,ist_p00)) / (amrex::Math::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +amrex::Math::abs(sten(i-1,j  ,k,ist_pp0)) + eps);
            Real w1p = amrex::Math::abs(sten(i  ,j,k,ist_p00)) / (amrex::Math::abs(sten(i  ,j-1,k,ist_pp0))
                                                         +amrex::Math::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real w2m = amrex::Math::abs(sten(i,j-1,k,ist_0p0)) / (amrex::Math::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +amrex::Math::abs(sten(i  ,j-1,k,ist_pp0)) + eps);
            Real w2p = amrex::Math::abs(sten(i,j  ,k,ist_0p0)) / (amrex::Math::abs(sten(i-1,j  ,k,ist_pp0))
                                                         +amrex::Math::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real wmm = amrex::Math::abs(sten(i-1,j-1,k,ist_pp0)) * (1.0 + w1m + w2m);
            Real wpm = amrex::Math::abs(sten(i  ,j-1,k,ist_pp0)) * (1.0 + w1p + w2m);
            Real wmp = amrex::Math::abs(sten(i-1,j  ,k,ist_pp0)) * (1.0 + w1m + w2p);
            Real wpp = amrex::Math::abs(sten(i  ,j  ,k,ist_pp0)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc+1,kc) + wpp*crse(ic+1,jc+1,kc))
                / (wmm+wpm+wmp+wpp+eps);
        } else {
            Real wmmm = 1.0;
            Real wpmm = 1.0;
            Real wmpm = 1.0;
            Real wppm = 1.0;
            Real wmmp = 1.0;
            Real wpmp = 1.0;
            Real wmpp = 1.0;
            Real wppp = 1.0;

            Real wtmp = amrex::Math::abs(sten(i-1,j,k,ist_p00)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j-1,k  ,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_p00)) /
                ( amrex::Math::abs(sten(i,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i,j  ,k-1,ist_ppp))
                + amrex::Math::abs(sten(i,j-1,k  ,ist_ppp))
                + amrex::Math::abs(sten(i,j  ,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j-1,k,ist_0p0)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j-1,k  ,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_0p0)) /
                ( amrex::Math::abs(sten(i-1,j,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j,k  ,ist_ppp))
                + amrex::Math::abs(sten(i  ,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k-1,ist_00p)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_00p)) /
                ( amrex::Math::abs(sten(i-1,j-1,k,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k,ist_ppp))
                + amrex::Math::abs(sten(i  ,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = amrex::Math::abs(sten(i-1,j-1,k,ist_pp0)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmmp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j-1,k,ist_pp0)) /
                ( amrex::Math::abs(sten(i,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i,j-1,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wpmp += wtmp;

            wtmp = amrex::Math::abs(sten(i-1,j,k,ist_pp0)) /
                ( amrex::Math::abs(sten(i-1,j,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wmpp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_pp0)) /
                ( amrex::Math::abs(sten(i,j,k-1,ist_ppp))
                + amrex::Math::abs(sten(i,j,k  ,ist_ppp)) + eps);
            wppm += wtmp;
            wppp += wtmp;

            wtmp = amrex::Math::abs(sten(i-1,j,k-1,ist_p0p)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k-1,ist_p0p)) /
                ( amrex::Math::abs(sten(i,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i,j  ,k-1,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;

            wtmp = amrex::Math::abs(sten(i-1,j,k,ist_p0p)) /
                ( amrex::Math::abs(sten(i-1,j-1,k,ist_ppp))
                + amrex::Math::abs(sten(i-1,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_p0p)) /
                ( amrex::Math::abs(sten(i,j-1,k,ist_ppp))
                + amrex::Math::abs(sten(i,j  ,k,ist_ppp)) + eps);
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j-1,k-1,ist_0pp)) /
                ( amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k-1,ist_0pp)) /
                ( amrex::Math::abs(sten(i-1,j,k-1,ist_ppp))
                + amrex::Math::abs(sten(i  ,j,k-1,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = amrex::Math::abs(sten(i,j-1,k,ist_0pp)) /
                ( amrex::Math::abs(sten(i-1,j-1,k,ist_ppp))
                + amrex::Math::abs(sten(i  ,j-1,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = amrex::Math::abs(sten(i,j,k,ist_0pp)) /
                ( amrex::Math::abs(sten(i-1,j,k,ist_ppp))
                + amrex::Math::abs(sten(i  ,j,k,ist_ppp)) + eps);
            wmpp += wtmp;
            wppp += wtmp;

            wmmm *= amrex::Math::abs(sten(i-1,j-1,k-1,ist_ppp));
            wpmm *= amrex::Math::abs(sten(i  ,j-1,k-1,ist_ppp));
            wmpm *= amrex::Math::abs(sten(i-1,j  ,k-1,ist_ppp));
            wppm *= amrex::Math::abs(sten(i  ,j  ,k-1,ist_ppp));
            wmmp *= amrex::Math::abs(sten(i-1,j-1,k  ,ist_ppp));
            wpmp *= amrex::Math::abs(sten(i  ,j-1,k  ,ist_ppp));
            wmpp *= amrex::Math::abs(sten(i-1,j  ,k  ,ist_ppp));
            wppp *= amrex::Math::abs(sten(i  ,j  ,k  ,ist_ppp));
            fv = (wmmm*crse(ic,jc  ,kc  ) + wpmm*crse(ic+1,jc  ,kc  )
                  + wmpm*crse(ic,jc+1,kc  ) + wppm*crse(ic+1,jc+1,kc  )
                  + wmmp*crse(ic,jc  ,kc+1) + wpmp*crse(ic+1,jc  ,kc+1)
                  + wmpp*crse(ic,jc+1,kc+1) + wppp*crse(ic+1,jc+1,kc+1))
                / (wmmm + wpmm + wmpm + wppm + wmmp + wpmp + wmpp + wppp + eps);
        }

        fine(i,j,k) += fv;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine, Array4<Real const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {

        Real cv = fine(ii,jj,kk);

        // ************************************
        // Adding fine(ii-1,jj,kk)
        // ************************************

        Real sten_lo = amrex::Math::abs(sten(ii-2,jj,kk,ist_p00));
        Real sten_hi = amrex::Math::abs(sten(ii-1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii-1,jj,kk);
        } else {
            cv += fine(ii-1,jj,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii+1,jj,kk)
        // ************************************

        sten_lo = amrex::Math::abs(sten(ii  ,jj,kk,ist_p00));
        sten_hi = amrex::Math::abs(sten(ii+1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii+1,jj,kk);
        } else {
            cv += fine(ii+1,jj,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj-1,kk)
        // ************************************

        sten_lo = amrex::Math::abs(sten(ii,jj-2,kk,ist_0p0));
        sten_hi = amrex::Math::abs(sten(ii,jj-1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj-1,kk);
        } else {
            cv += fine(ii,jj-1,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj+1,kk)
        // ************************************

        sten_lo = amrex::Math::abs(sten(ii,jj  ,kk,ist_0p0));
        sten_hi = amrex::Math::abs(sten(ii,jj+1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj+1,kk);
        } else {
            cv += fine(ii,jj+1,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk-1)
        // ************************************

        sten_lo = amrex::Math::abs(sten(ii,jj,kk-2,ist_00p));
        sten_hi = amrex::Math::abs(sten(ii,jj,kk-1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk-1);
        } else {
            cv += fine(ii,jj,kk-1)*sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk+1)
        // ************************************

        sten_lo = amrex::Math::abs(sten(ii,jj,kk  ,ist_00p));
        sten_hi = amrex::Math::abs(sten(ii,jj,kk+1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk+1);
        } else {
            cv += fine(ii,jj,kk+1)*sten_lo  / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii-1,jj-1,kk)
        // ************************************

        // keven
        Real w1m = amrex::Math::abs(sten(ii-2,jj-1,kk,ist_p00))
            / (    amrex::Math::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +amrex::Math::abs(sten(ii-2,jj-1,kk,ist_pp0)) + eps);
        Real w1p = amrex::Math::abs(sten(ii-1,jj-1,kk,ist_p00))
            / (    amrex::Math::abs(sten(ii-1,jj-2,kk,ist_pp0))
                  +amrex::Math::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real w2m = amrex::Math::abs(sten(ii-1,jj-2,kk,ist_0p0))
            / (    amrex::Math::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +amrex::Math::abs(sten(ii-1,jj-2,kk,ist_pp0)) + eps);
        Real w2p = amrex::Math::abs(sten(ii-1,jj-1,kk,ist_0p0))
            / (    amrex::Math::abs(sten(ii-2,jj-1,kk,ist_pp0))
                  +amrex::Math::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real wmm = amrex::Math::abs(sten(ii-2,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        Real wpm = amrex::Math::abs(sten(ii-1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        Real wmp = amrex::Math::abs(sten(ii-2,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        Real wpp = amrex::Math::abs(sten(ii-1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj-1,kk)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj-1,kk)
        // ************************************

        w1m = amrex::Math::abs(sten(ii  ,jj-1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +amrex::Math::abs(sten(ii  ,jj-1,kk,ist_pp0)) + eps);
        w1p = amrex::Math::abs(sten(ii+1,jj-1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii+1,jj-2,kk,ist_pp0))
             +amrex::Math::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        w2m = amrex::Math::abs(sten(ii+1,jj-2,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +amrex::Math::abs(sten(ii+1,jj-2,kk,ist_pp0)) + eps);
        w2p = amrex::Math::abs(sten(ii+1,jj-1,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii  ,jj-1,kk,ist_pp0))
             +amrex::Math::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        wmm = amrex::Math::abs(sten(ii  ,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii+1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii  ,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii+1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj-1,kk)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj+1,kk)
        // ************************************

        w1m = amrex::Math::abs(sten(ii-2,jj+1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +amrex::Math::abs(sten(ii-2,jj+1,kk,ist_pp0)) + eps);
        w1p = amrex::Math::abs(sten(ii-1,jj+1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii-1,jj  ,kk,ist_pp0))
             +amrex::Math::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        w2m = amrex::Math::abs(sten(ii-1,jj  ,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +amrex::Math::abs(sten(ii-1,jj  ,kk,ist_pp0)) + eps);
        w2p = amrex::Math::abs(sten(ii-1,jj+1,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii-2,jj+1,kk,ist_pp0))
             +amrex::Math::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        wmm = amrex::Math::abs(sten(ii-2,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii-1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii-2,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii-1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj+1,kk)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj+1,kk)
        // ************************************

        w1m = amrex::Math::abs(sten(ii  ,jj+1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii  ,jj+1,kk,ist_pp0))
             +amrex::Math::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w1p = amrex::Math::abs(sten(ii+1,jj+1,kk,ist_p00))
           / (amrex::Math::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +amrex::Math::abs(sten(ii+1,jj  ,kk,ist_pp0)) + eps);
        w2m = amrex::Math::abs(sten(ii+1,jj  ,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii+1,jj  ,kk,ist_pp0))
             +amrex::Math::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w2p = amrex::Math::abs(sten(ii+1,jj+1,kk,ist_0p0))
           / (amrex::Math::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +amrex::Math::abs(sten(ii  ,jj+1,kk,ist_pp0)) + eps);
        wmm = amrex::Math::abs(sten(ii  ,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii+1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii  ,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii+1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj+1,kk)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk-1)
        // ************************************

        // jeven
        w1m = amrex::Math::abs(sten(ii-2,jj,kk-1,ist_p00))
           / (amrex::Math::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii-2,jj,kk-1,ist_p0p)) + eps);
        w1p = amrex::Math::abs(sten(ii-1,jj,kk-1,ist_p00))
           / (amrex::Math::abs(sten(ii-1,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        w2m = amrex::Math::abs(sten(ii-1,jj,kk-2,ist_00p))
           / (amrex::Math::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk-2,ist_p0p)) + eps);
        w2p = amrex::Math::abs(sten(ii-1,jj,kk-1,ist_00p))
           / (amrex::Math::abs(sten(ii-2,jj,kk-1,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        wmm = amrex::Math::abs(sten(ii-2,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii-1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii-2,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii-1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk-1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii  ,jj,kk-1,ist_p00))
           / (amrex::Math::abs(sten(ii  ,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        w1p = amrex::Math::abs(sten(ii+1,jj,kk-1,ist_p00))
           / (amrex::Math::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii+1,jj,kk-1,ist_p0p)) + eps);
        w2m = amrex::Math::abs(sten(ii+1,jj,kk-2,ist_00p))
           / (amrex::Math::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk-2,ist_p0p)) + eps);
        w2p = amrex::Math::abs(sten(ii+1,jj,kk-1,ist_00p))
           / (amrex::Math::abs(sten(ii+1,jj,kk-1,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        wmm = amrex::Math::abs(sten(ii  ,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii+1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii  ,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii+1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk+1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii-2,jj,kk+1,ist_p00))
           / (amrex::Math::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +amrex::Math::abs(sten(ii-2,jj,kk  ,ist_p0p)) + eps);
        w1p = amrex::Math::abs(sten(ii-1,jj,kk+1,ist_p00))
           / (amrex::Math::abs(sten(ii-1,jj,kk+1,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2m = amrex::Math::abs(sten(ii-1,jj,kk  ,ist_00p))
           / (amrex::Math::abs(sten(ii-2,jj,kk  ,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2p = amrex::Math::abs(sten(ii-1,jj,kk+1,ist_00p))
           / (amrex::Math::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +amrex::Math::abs(sten(ii-1,jj,kk+1,ist_p0p)) + eps);
        wmm = amrex::Math::abs(sten(ii-2,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii-1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii-2,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii-1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk+1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii  ,jj,kk+1,ist_p00))
           / (amrex::Math::abs(sten(ii  ,jj,kk+1,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w1p = amrex::Math::abs(sten(ii+1,jj,kk+1,ist_p00))
           / (amrex::Math::abs(sten(ii+1,jj,kk+1,ist_p0p))
              +amrex::Math::abs(sten(ii+1,jj,kk  ,ist_p0p)) + eps);
        w2m = amrex::Math::abs(sten(ii+1,jj,kk  ,ist_00p))
           / (amrex::Math::abs(sten(ii+1,jj,kk  ,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w2p = amrex::Math::abs(sten(ii+1,jj,kk+1,ist_00p))
           / (amrex::Math::abs(sten(ii+1,jj,kk+1,ist_p0p))
             +amrex::Math::abs(sten(ii  ,jj,kk+1,ist_p0p)) + eps);
        wmm = amrex::Math::abs(sten(ii  ,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii+1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii  ,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii+1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk-1)
        // ************************************

        // ieven
        w1m = amrex::Math::abs(sten(ii,jj-2,kk-1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-2,kk-1,ist_0pp)) + eps);
        w2m = amrex::Math::abs(sten(ii,jj-1,kk-2,ist_00p))
           / (amrex::Math::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk-2,ist_0pp)) + eps);
        w1p = amrex::Math::abs(sten(ii,jj-1,kk-1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj-1,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        w2p = amrex::Math::abs(sten(ii,jj-1,kk-1,ist_00p))
           / (amrex::Math::abs(sten(ii,jj-2,kk-1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        wmm = amrex::Math::abs(sten(ii,jj-2,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii,jj-1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii,jj-2,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii,jj-1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk-1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii,jj  ,kk-1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj  ,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        w1p = amrex::Math::abs(sten(ii,jj+1,kk-1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj+1,kk-1,ist_0pp)) + eps);
        w2m = amrex::Math::abs(sten(ii,jj+1,kk-2,ist_00p))
           / (amrex::Math::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk-2,ist_0pp)) + eps);
        w2p = amrex::Math::abs(sten(ii,jj+1,kk-1,ist_00p))
           / (amrex::Math::abs(sten(ii,jj+1,kk-1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        wmm = amrex::Math::abs(sten(ii,jj  ,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii,jj+1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii,jj  ,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii,jj+1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk+1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii,jj-2,kk+1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-2,kk  ,ist_0pp)) + eps);
        w1p = amrex::Math::abs(sten(ii,jj-1,kk+1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj-1,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2m = amrex::Math::abs(sten(ii,jj-1,kk  ,ist_00p))
           / (amrex::Math::abs(sten(ii,jj-2,kk  ,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2p = amrex::Math::abs(sten(ii,jj-1,kk+1,ist_00p))
           / (amrex::Math::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj-1,kk+1,ist_0pp)) + eps);
        wmm = amrex::Math::abs(sten(ii,jj-2,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii,jj-1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii,jj-2,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii,jj-1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk+1)
        // ************************************

        w1m = amrex::Math::abs(sten(ii,jj  ,kk+1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj  ,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w1p = amrex::Math::abs(sten(ii,jj+1,kk+1,ist_0p0))
           / (amrex::Math::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj+1,kk  ,ist_0pp)) + eps);
        w2m = amrex::Math::abs(sten(ii,jj+1,kk  ,ist_00p))
           / (amrex::Math::abs(sten(ii,jj+1,kk  ,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w2p = amrex::Math::abs(sten(ii,jj+1,kk+1,ist_00p))
           / (amrex::Math::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +amrex::Math::abs(sten(ii,jj  ,kk+1,ist_0pp)) + eps);
        wmm = amrex::Math::abs(sten(ii,jj  ,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = amrex::Math::abs(sten(ii,jj+1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = amrex::Math::abs(sten(ii,jj  ,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = amrex::Math::abs(sten(ii,jj+1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine at corners
        // ************************************

        Real wmmm = 1.0
            + amrex::Math::abs(sten(ii  ,jj+1,kk+1,ist_p00)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj  ,kk+1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj+1,kk  ,ist_00p)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj+1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj  ,kk+1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj+1,kk  ,ist_p0p)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj  ,kk  ,ist_0pp)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk  ,ist_ppp)) + eps);
        wmmm *= amrex::Math::abs(sten(ii,jj,kk,ist_ppp));
        cv += wmmm*fine(ii+1,jj+1,kk+1)*sten(ii+1,jj+1,kk+1,ist_inv);

        Real wpmm = 1.0
            + amrex::Math::abs(sten(ii-1,jj+1,kk+1,ist_p00)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk+1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj  ,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj+1,kk  ,ist_00p)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj+1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk+1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj+1,kk  ,ist_p0p)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_0pp)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk  ,ist_ppp)) + eps);
        wpmm *= amrex::Math::abs(sten(ii-1,jj,kk,ist_ppp));
        cv += wpmm*fine(ii-1,jj+1,kk+1)*sten(ii-1,jj+1,kk+1,ist_inv);

        Real wmpm = 1.0
            + amrex::Math::abs(sten(ii  ,jj-1,kk+1,ist_p00)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-2,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk+1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk  ,ist_00p)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj-1,kk+1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_p0p)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk  ,ist_0pp)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps);
        wmpm *= amrex::Math::abs(sten(ii  ,jj-1,kk  ,ist_ppp));
        cv += wmpm*fine(ii+1,jj-1,kk+1)*sten(ii+1,jj-1,kk+1,ist_inv);

        Real wppm = 1.0
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_p00)) /
            ( amrex::Math::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-2,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj-1,kk+1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_00p)) /
            ( amrex::Math::abs(sten(ii-2,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_p0p)) /
            ( amrex::Math::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_0pp)) /
            ( amrex::Math::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps);
        wppm *= amrex::Math::abs(sten(ii-1,jj-1,kk  ,ist_ppp));
        cv += wppm*fine(ii-1,jj-1,kk+1)*sten(ii-1,jj-1,kk+1,ist_inv);

        Real wmmp = 1.0
            + amrex::Math::abs(sten(ii  ,jj+1,kk-1,ist_p00)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj  ,kk-1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj+1,kk-1,ist_00p)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj+1,kk-1,ist_p0p)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj  ,kk-1,ist_0pp)) /
            ( amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps);
        wmmp *= amrex::Math::abs(sten(ii  ,jj  ,kk-1,ist_ppp));
        cv += wmmp*fine(ii+1,jj+1,kk-1)*sten(ii+1,jj+1,kk-1,ist_inv);

        Real wpmp = 1.0
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_p00)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_00p)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj+1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_p0p)) /
            ( amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_0pp)) /
            ( amrex::Math::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps);
        wpmp *= amrex::Math::abs(sten(ii-1,jj  ,kk-1,ist_ppp));
        cv += wpmp*fine(ii-1,jj+1,kk-1)*sten(ii-1,jj+1,kk-1,ist_inv);

        Real wmpp = 1.0
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_p00)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_00p)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_p0p)) /
            ( amrex::Math::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_0pp)) /
            ( amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps);
        wmpp *= amrex::Math::abs(sten(ii  ,jj-1,kk-1,ist_ppp));
        cv += wmpp*fine(ii+1,jj-1,kk-1)*sten(ii+1,jj-1,kk-1,ist_inv);

        Real wppp = 1.0
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_p00)) /
            ( amrex::Math::abs(sten(ii-1,jj-2,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_0p0)) /
            ( amrex::Math::abs(sten(ii-2,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_00p)) /
            ( amrex::Math::abs(sten(ii-2,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_pp0)) /
            ( amrex::Math::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_p0p)) /
            ( amrex::Math::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_0pp)) /
            ( amrex::Math::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps);
        wppp *= amrex::Math::abs(sten(ii-1,jj-1,kk-1,ist_ppp));
        cv += wppp*fine(ii-1,jj-1,kk-1)*sten(ii-1,jj-1,kk-1,ist_inv);

        crse(i,j,k) = cv * 0.125;
    }
}

#ifdef AMREX_USE_EB

namespace {

    constexpr int i_S_x     = 0;
    constexpr int i_S_y     = 1;
    constexpr int i_S_z     = 2;
    constexpr int i_S_x2    = 3;
    constexpr int i_S_y2    = 4;
    constexpr int i_S_z2    = 5;
    constexpr int i_S_x_y   = 6;
    constexpr int i_S_x_z   = 7;
    constexpr int i_S_y_z   = 8;
    constexpr int i_S_x2_y  = 9;
    constexpr int i_S_x2_z  = 10;
    constexpr int i_S_x_y2  = 11;
    constexpr int i_S_y2_z  = 12;
    constexpr int i_S_x_z2  = 13;
    constexpr int i_S_y_z2  = 14;
    constexpr int i_S_x2_y2 = 15;
    constexpr int i_S_x2_z2 = 16;
    constexpr int i_S_y2_z2 = 17;
    constexpr int i_S_xyz   = 18;
    constexpr int n_Sintg   = 19;

    constexpr int i_c_xmym = 0;
    constexpr int i_c_xmyb = 1;
    constexpr int i_c_xmyp = 2;
    constexpr int i_c_xbym = 3;
    constexpr int i_c_xbyb = 4;
    constexpr int i_c_xbyp = 5;
    constexpr int i_c_xpym = 6;
    constexpr int i_c_xpyb = 7;
    constexpr int i_c_xpyp = 8;
    constexpr int i_c_xmzm = 9;
    constexpr int i_c_xmzb = 10;
    constexpr int i_c_xmzp = 11;
    constexpr int i_c_xbzm = 12;
    constexpr int i_c_xbzb = 13;
    constexpr int i_c_xbzp = 14;
    constexpr int i_c_xpzm = 15;
    constexpr int i_c_xpzb = 16;
    constexpr int i_c_xpzp = 17;
    constexpr int i_c_ymzm = 18;
    constexpr int i_c_ymzb = 19;
    constexpr int i_c_ymzp = 20;
    constexpr int i_c_ybzm = 21;
    constexpr int i_c_ybzb = 22;
    constexpr int i_c_ybzp = 23;
    constexpr int i_c_ypzm = 24;
    constexpr int i_c_ypzb = 25;
    constexpr int i_c_ypzp = 26;
    constexpr int n_conn = 27;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_connection (int i, int j, int k, Array4<Real> const& conn,
                             Array4<Real const> const& intg, Array4<Real const> const& vol,
                             Array4<EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,k).isCovered()) {
        for (int n = 0; n < n_conn; ++n) conn(i,j,k,n) = Real(0.);
    } else if (flag(i,j,k).isRegular() or vol(i,j,k) >= almostone) {
        for (int n = 0; n < n_conn; ++n) conn(i,j,k,n) = Real(1.);
    } else {
        // Scaled by 9
        conn(i,j,k,i_c_xmym) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(9.)*( intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xmyb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + Real(18.)*( intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmyp) =  Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(9.)*(-intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xbym) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(18.)*(intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 36
        conn(i,j,k,i_c_xbyb) = Real(2.25)*vol(i,j,k)
            + Real(9.)*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + Real(36.)*intg(i,j,k,i_S_x2_y2);

        // Scaled by 18
        conn(i,j,k,i_c_xbyp) =  Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(18.)*(-intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpym) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(9.)*(-intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xpyb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + Real(18.)*(-intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpyp) =  Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + Real(9.)*( intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                      +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzm) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_x) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(intg(i,j,k,i_S_x_z) - intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xmzb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + Real(18.)*(intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzp) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_x  ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(-intg(i,j,k,i_S_x_z  ) + intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzm) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(18.)*(intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzb) = Real(2.25)*vol(i,j,k)
            + Real(9.)*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + Real(36.)*intg(i,j,k,i_S_x2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_xbzp) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(18.)*(-intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzm) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(-intg(i,j,k,i_S_x_z ) - intg(i,j,k,i_S_x2_z)
                     +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xpzb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_x   ) + intg(i,j,k,i_S_x2   ) - intg(i,j,k,i_S_z2))
            + Real(18.)*(-intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzp) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + Real(9.)*( intg(i,j,k,i_S_x_z ) + intg(i,j,k,i_S_x2_z)
                      +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzm) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(intg(i,j,k,i_S_y_z) - intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ymzb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_y) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + Real(18.)*(intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzp) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*(-intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(-intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ybzm) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(18.)*(intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 36
        conn(i,j,k,i_c_ybzb) = Real(2.25)*vol(i,j,k)
            + Real(9.)*(-intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
                     + Real(36.)*intg(i,j,k,i_S_y2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_ybzp) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(18.)*(-intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzm) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_y ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(9.)*(-intg(i,j,k,i_S_y_z ) - intg(i,j,k,i_S_y2_z)
                     +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ypzb) = Real(1.125)*vol(i,j,k)
            + Real(4.5)*( intg(i,j,k,i_S_y   ) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + Real(18.)*(-intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzp) = Real(0.5625)*vol(i,j,k)
            + Real(2.25)*( intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + Real(9.)*( intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                      +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_eb (int i, int j, int k, Array4<Real> const& sten,
                             Array4<Real const> const& sig, Array4<Real const> const& conn,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = Real(1./36.)*dxinv[0]*dxinv[0];
    Real facy = Real(1./36.)*dxinv[1]*dxinv[1];
    Real facz = Real(1./36.)*dxinv[2]*dxinv[2];

    // i+1,j,k
    sten(i,j,k,ist_p00) = (
        sig(i,j  ,k  )*(Real(4.)*facx*conn(i,j  ,k  ,i_c_ymzm) - Real(2.)*facy*conn(i,j  ,k  ,i_c_xbzm) - Real(2.)*facz*conn(i,j  ,k  ,i_c_xbym) ) +
        sig(i,j-1,k  )*(Real(4.)*facx*conn(i,j-1,k  ,i_c_ypzm) - Real(2.)*facy*conn(i,j-1,k  ,i_c_xbzm) - Real(2.)*facz*conn(i,j-1,k  ,i_c_xbyp) ) +
        sig(i,j  ,k-1)*(Real(4.)*facx*conn(i,j  ,k-1,i_c_ymzp) - Real(2.)*facy*conn(i,j  ,k-1,i_c_xbzp) - Real(2.)*facz*conn(i,j  ,k-1,i_c_xbym) ) +
        sig(i,j-1,k-1)*(Real(4.)*facx*conn(i,j-1,k-1,i_c_ypzp) - Real(2.)*facy*conn(i,j-1,k-1,i_c_xbzp) - Real(2.)*facz*conn(i,j-1,k-1,i_c_xbyp) ) );

    // i,j+1,k
    sten(i,j,k,ist_0p0) = (
        sig(i  ,j,k  )*(Real(-2.)*facx*conn(i  ,j,k  ,i_c_ybzm) + Real(4.)*facy*conn(i  ,j,k  ,i_c_xmzm) - Real(2.)*facz*conn(i  ,j,k  ,i_c_xmyb) ) +
        sig(i-1,j,k  )*(Real(-2.)*facx*conn(i-1,j,k  ,i_c_ybzm) + Real(4.)*facy*conn(i-1,j,k  ,i_c_xpzm) - Real(2.)*facz*conn(i-1,j,k  ,i_c_xpyb) ) +
        sig(i  ,j,k-1)*(Real(-2.)*facx*conn(i  ,j,k-1,i_c_ybzp) + Real(4.)*facy*conn(i  ,j,k-1,i_c_xmzp) - Real(2.)*facz*conn(i  ,j,k-1,i_c_xmyb) ) +
        sig(i-1,j,k-1)*(Real(-2.)*facx*conn(i-1,j,k-1,i_c_ybzp) + Real(4.)*facy*conn(i-1,j,k-1,i_c_xpzp) - Real(2.)*facz*conn(i-1,j,k-1,i_c_xpyb) ) );

    // i,j,k+1
    sten(i,j,k,ist_00p) = (
        sig(i  ,j  ,k)*(Real(-2.)*facx*conn(i  ,j  ,k,i_c_ymzb) - Real(2.)*facy*conn(i  ,j  ,k,i_c_xmzb) + Real(4.)*facz*conn(i  ,j  ,k,i_c_xmym) ) +
        sig(i-1,j  ,k)*(Real(-2.)*facx*conn(i-1,j  ,k,i_c_ymzb) - Real(2.)*facy*conn(i-1,j  ,k,i_c_xpzb) + Real(4.)*facz*conn(i-1,j  ,k,i_c_xpym) ) +
        sig(i  ,j-1,k)*(Real(-2.)*facx*conn(i  ,j-1,k,i_c_ypzb) - Real(2.)*facy*conn(i  ,j-1,k,i_c_xmzb) + Real(4.)*facz*conn(i  ,j-1,k,i_c_xmyp) ) +
        sig(i-1,j-1,k)*(Real(-2.)*facx*conn(i-1,j-1,k,i_c_ypzb) - Real(2.)*facy*conn(i-1,j-1,k,i_c_xpzb) + Real(4.)*facz*conn(i-1,j-1,k,i_c_xpyp) ) );

    // i+1,j+1,k
    sten(i,j,k,ist_pp0) = (
        sig(i,j,k  )*(Real(2.)*facx*conn(i,j,k  ,i_c_ybzm) + Real(2.)*facy*conn(i,j,k  ,i_c_xbzm) - facz*conn(i,j,k  ,i_c_xbyb) ) +
        sig(i,j,k-1)*(Real(2.)*facx*conn(i,j,k-1,i_c_ybzp) + Real(2.)*facy*conn(i,j,k-1,i_c_xbzp) - facz*conn(i,j,k-1,i_c_xbyb) ) );

    // i+1,j,k+1
    sten(i,j,k,ist_p0p) = (
        sig(i,j,k  )*(Real(2.)*facx*conn(i,j,k  ,i_c_ymzb) - facy*conn(i,j,k  ,i_c_xbzb) + Real(2.)*facz*conn(i,j,k  ,i_c_xbym) ) +
        sig(i,j-1,k)*(Real(2.)*facx*conn(i,j-1,k,i_c_ypzb) - facy*conn(i,j-1,k,i_c_xbzb) + Real(2.)*facz*conn(i,j-1,k,i_c_xbyp) ) );

    // i,j+1,k+1
    sten(i,j,k,ist_0pp) = (
        sig(i  ,j,k)*(-facx*conn(i  ,j,k,i_c_ybzb) + Real(2.)*facy*conn(i  ,j,k,i_c_xmzb) + Real(2.)*facz*conn(i  ,j,k,i_c_xmyb) ) +
        sig(i-1,j,k)*(-facx*conn(i-1,j,k,i_c_ybzb) + Real(2.)*facy*conn(i-1,j,k,i_c_xpzb) + Real(2.)*facz*conn(i-1,j,k,i_c_xpyb) ) );

    // i+1,j+1,k+1
    sten(i,j,k,ist_ppp) = sig(i,j,k) * (facx*conn(i,j,k,i_c_ybzb) + facy*conn(i,j,k,i_c_xbzb) + facz*conn(i,j,k,i_c_xbyb) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_eb (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Box const& nodal_domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    Real facx = Real(0.25)*dxinv[0];
    Real facy = Real(0.25)*dxinv[1];
    Real facz = Real(0.25)*dxinv[2];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (!msk(i,j,k)) {

        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;
        Real zero_klo = 1.0;
        Real zero_khi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0;
        if ((bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)
            and k == domlo.z) zero_klo = 0.0;
        if ((bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)
            and k == domhi.z) zero_khi = 0.0;

        rhs(i,j,k) = facx*(
            vel(i-1,j-1,k  ,0)*(    -vfrac(i-1,j-1,k  )
                               -Real(2.)*intg(i-1,j-1,k  ,i_S_y)
                               +Real(2.)*intg(i-1,j-1,k  ,i_S_z)
                               +Real(4.)*intg(i-1,j-1,k  ,i_S_y_z))*zero_jlo*zero_khi
           +vel(i  ,j-1,k  ,0)*(     vfrac(i  ,j-1,k  )
                               +Real(2.)*intg(i  ,j-1,k  ,i_S_y)
                               -Real(2.)*intg(i  ,j-1,k  ,i_S_z)
                               -Real(4.)*intg(i  ,j-1,k  ,i_S_y_z))*zero_jlo*zero_khi
           +vel(i-1,j  ,k  ,0)*(    -vfrac(i-1,j  ,k  )
                               +Real(2.)*intg(i-1,j  ,k  ,i_S_y)
                               +Real(2.)*intg(i-1,j  ,k  ,i_S_z)
                               -Real(4.)*intg(i-1,j  ,k  ,i_S_y_z))*zero_jhi*zero_khi
           +vel(i  ,j  ,k  ,0)*(     vfrac(i  ,j  ,k  )
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_y)
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_z)
                               +Real(4.)*intg(i  ,j  ,k  ,i_S_y_z))*zero_jhi*zero_khi
           +vel(i-1,j-1,k-1,0)*(    -vfrac(i-1,j-1,k-1)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_y)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_z)
                               -Real(4.)*intg(i-1,j-1,k-1,i_S_y_z))*zero_jlo*zero_klo
           +vel(i  ,j-1,k-1,0)*(     vfrac(i  ,j-1,k-1)
                               +Real(2.)*intg(i  ,j-1,k-1,i_S_y)
                               +Real(2.)*intg(i  ,j-1,k-1,i_S_z)
                               +Real(4.)*intg(i  ,j-1,k-1,i_S_y_z))*zero_jlo*zero_klo
           +vel(i-1,j  ,k-1,0)*(    -vfrac(i-1,j  ,k-1)
                               +Real(2.)*intg(i-1,j  ,k-1,i_S_y)
                               -Real(2.)*intg(i-1,j  ,k-1,i_S_z)
                               +Real(4.)*intg(i-1,j  ,k-1,i_S_y_z))*zero_jhi*zero_klo
           +vel(i  ,j  ,k-1,0)*(     vfrac(i  ,j  ,k-1)
                               -Real(2.)*intg(i  ,j  ,k-1,i_S_y)
                               +Real(2.)*intg(i  ,j  ,k-1,i_S_z)
                               -Real(4.)*intg(i  ,j  ,k-1,i_S_y_z))*zero_jhi*zero_klo )
            + facy*(
            vel(i-1,j-1,k  ,1)*(    -vfrac(i-1,j-1,k  )
                               -Real(2.)*intg(i-1,j-1,k  ,i_S_x)
                               +Real(2.)*intg(i-1,j-1,k  ,i_S_z)
                               +Real(4.)*intg(i-1,j-1,k  ,i_S_x_z))*zero_ilo*zero_khi
           +vel(i  ,j-1,k  ,1)*(    -vfrac(i  ,j-1,k  )
                               +Real(2.)*intg(i  ,j-1,k  ,i_S_x)
                               +Real(2.)*intg(i  ,j-1,k  ,i_S_z)
                               -Real(4.)*intg(i  ,j-1,k  ,i_S_x_z))*zero_ihi*zero_khi
           +vel(i-1,j  ,k  ,1)*(     vfrac(i-1,j  ,k  )
                               +Real(2.)*intg(i-1,j  ,k  ,i_S_x)
                               -Real(2.)*intg(i-1,j  ,k  ,i_S_z)
                               -Real(4.)*intg(i-1,j  ,k  ,i_S_x_z))*zero_ilo*zero_khi
           +vel(i  ,j  ,k  ,1)*(     vfrac(i  ,j  ,k  )
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_x)
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_z)
                               +Real(4.)*intg(i  ,j  ,k  ,i_S_x_z))*zero_ihi*zero_khi
           +vel(i-1,j-1,k-1,1)*(    -vfrac(i-1,j-1,k-1)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_x)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_z)
                               -Real(4.)*intg(i-1,j-1,k-1,i_S_x_z))*zero_ilo*zero_klo
           +vel(i  ,j-1,k-1,1)*(    -vfrac(i  ,j-1,k-1)
                               +Real(2.)*intg(i  ,j-1,k-1,i_S_x)
                               -Real(2.)*intg(i  ,j-1,k-1,i_S_z)
                               +Real(4.)*intg(i  ,j-1,k-1,i_S_x_z))*zero_ihi*zero_klo
           +vel(i-1,j  ,k-1,1)*(     vfrac(i-1,j  ,k-1)
                               +Real(2.)*intg(i-1,j  ,k-1,i_S_x)
                               +Real(2.)*intg(i-1,j  ,k-1,i_S_z)
                               +Real(4.)*intg(i-1,j  ,k-1,i_S_x_z))*zero_ilo*zero_klo
           +vel(i  ,j  ,k-1,1)*(     vfrac(i  ,j  ,k-1)
                               -Real(2.)*intg(i  ,j  ,k-1,i_S_x)
                               +Real(2.)*intg(i  ,j  ,k-1,i_S_z)
                               -Real(4.)*intg(i  ,j  ,k-1,i_S_x_z))*zero_ihi*zero_klo )
            + facz*(
            vel(i-1,j-1,k  ,2)*(     vfrac(i-1,j-1,k  )
                               +Real(2.)*intg(i-1,j-1,k  ,i_S_x)
                               +Real(2.)*intg(i-1,j-1,k  ,i_S_y)
                               +Real(4.)*intg(i-1,j-1,k  ,i_S_x_y))*zero_ilo*zero_jlo
           +vel(i  ,j-1,k  ,2)*(     vfrac(i  ,j-1,k  )
                               -Real(2.)*intg(i  ,j-1,k  ,i_S_x)
                               +Real(2.)*intg(i  ,j-1,k  ,i_S_y)
                               -Real(4.)*intg(i  ,j-1,k  ,i_S_x_y))*zero_ihi*zero_jlo
           +vel(i-1,j  ,k  ,2)*(     vfrac(i-1,j  ,k  )
                               +Real(2.)*intg(i-1,j  ,k  ,i_S_x)
                               -Real(2.)*intg(i-1,j  ,k  ,i_S_y)
                               -Real(4.)*intg(i-1,j  ,k  ,i_S_x_y))*zero_ilo*zero_jhi
           +vel(i  ,j  ,k  ,2)*(     vfrac(i  ,j  ,k  )
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_x)
                               -Real(2.)*intg(i  ,j  ,k  ,i_S_y)
                               +Real(4.)*intg(i  ,j  ,k  ,i_S_x_y))*zero_ihi*zero_jhi
           +vel(i-1,j-1,k-1,2)*(    -vfrac(i-1,j-1,k-1)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_x)
                               -Real(2.)*intg(i-1,j-1,k-1,i_S_y)
                               -Real(4.)*intg(i-1,j-1,k-1,i_S_x_y))*zero_ilo*zero_jlo
           +vel(i  ,j-1,k-1,2)*(    -vfrac(i  ,j-1,k-1)
                               +Real(2.)*intg(i  ,j-1,k-1,i_S_x)
                               -Real(2.)*intg(i  ,j-1,k-1,i_S_y)
                               +Real(4.)*intg(i  ,j-1,k-1,i_S_x_y))*zero_ihi*zero_jlo
           +vel(i-1,j  ,k-1,2)*(    -vfrac(i-1,j  ,k-1)
                               -Real(2.)*intg(i-1,j  ,k-1,i_S_x)
                               +Real(2.)*intg(i-1,j  ,k-1,i_S_y)
                               +Real(4.)*intg(i-1,j  ,k-1,i_S_x_y))*zero_ilo*zero_jhi
           +vel(i  ,j  ,k-1,2)*(    -vfrac(i  ,j  ,k-1)
                               +Real(2.)*intg(i  ,j  ,k-1,i_S_x)
                               +Real(2.)*intg(i  ,j  ,k-1,i_S_y)
                               -Real(4.)*intg(i  ,j  ,k-1,i_S_x_y))*zero_ihi*zero_jhi );
    } else {
        rhs(i,j,k) = Real(0.);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                        Array4<Real const> const& sig, Array4<Real const> const& vfrac,
                        Array4<Real const> const& intg, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (vfrac(i,j,k) == Real(0.)) {
        u(i,j,k,0) = u(i,j,k,1) = u(i,j,k,2) = Real(0.);
    } else {
        Real dpdx = Real(0.25)*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
        Real dpdy = Real(0.25)*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
        Real dpdz = Real(0.25)*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
                                +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));

        Real dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k);

        u(i,j,k,0) -= sig(i,j,k)*dxinv[0]*(dpdx + Real(0.5)*intg(i,j,k,i_S_y  )*dpp_xy +
                                                  Real(0.5)*intg(i,j,k,i_S_z  )*dpp_xz +
                                                            intg(i,j,k,i_S_y_z)*dpp_xyz );
        u(i,j,k,1) -= sig(i,j,k)*dxinv[1]*(dpdy + Real(0.5)*intg(i,j,k,i_S_x  )*dpp_xy +
                                                  Real(0.5)*intg(i,j,k,i_S_z  )*dpp_yz +
                                                            intg(i,j,k,i_S_x_z)*dpp_xyz );
        u(i,j,k,2) -= sig(i,j,k)*dxinv[2]*(dpdz + Real(0.5)*intg(i,j,k,i_S_x  )*dpp_xz +
                                                  Real(0.5)*intg(i,j,k,i_S_y  )*dpp_yz +
                                                            intg(i,j,k,i_S_x_y)*dpp_xyz );
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc_eb (int i, int j, int k, Array4<Real const> const& rhcc,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        return
                          rhcc(i  ,j  ,k  ) *
         ( Real(0.125) * vfrac(i  ,j  ,k  )
         + Real(0.25) * (-intg(i  ,j  ,k  ,i_S_x)
                         -intg(i  ,j  ,k  ,i_S_y)
                         -intg(i  ,j  ,k  ,i_S_z))
         + Real(0.5) * (  intg(i  ,j  ,k  ,i_S_x_y)
                         +intg(i  ,j  ,k  ,i_S_x_z)
                         +intg(i  ,j  ,k  ,i_S_y_z))
            +          ( -intg(i  ,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k  ) *
         ( Real(0.125) * vfrac(i-1,j  ,k  )
         + Real(0.25) * ( intg(i-1,j  ,k  ,i_S_x)
                         -intg(i-1,j  ,k  ,i_S_y)
                         -intg(i-1,j  ,k  ,i_S_z))
         + Real(0.5) * ( -intg(i-1,j  ,k  ,i_S_x_y)
                         -intg(i-1,j  ,k  ,i_S_x_z)
                         +intg(i-1,j  ,k  ,i_S_y_z))
            +          (  intg(i-1,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k  ) *
         ( Real(0.125) * vfrac(i  ,j-1,k  )
         + Real(0.25) * (-intg(i  ,j-1,k  ,i_S_x)
                         +intg(i  ,j-1,k  ,i_S_y)
                         -intg(i  ,j-1,k  ,i_S_z))
         + Real(0.5) * ( -intg(i  ,j-1,k  ,i_S_x_y)
                         +intg(i  ,j-1,k  ,i_S_x_z)
                         -intg(i  ,j-1,k  ,i_S_y_z))
            +          (  intg(i  ,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k  ) *
         ( Real(0.125) * vfrac(i-1,j-1,k  )
         + Real(0.25) * ( intg(i-1,j-1,k  ,i_S_x)
                         +intg(i-1,j-1,k  ,i_S_y)
                         -intg(i-1,j-1,k  ,i_S_z))
         + Real(0.5) * (  intg(i-1,j-1,k  ,i_S_x_y)
                         -intg(i-1,j-1,k  ,i_S_x_z)
                         -intg(i-1,j-1,k  ,i_S_y_z))
            +          ( -intg(i-1,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j  ,k-1) *
         ( Real(0.125) * vfrac(i  ,j  ,k-1)
         + Real(0.25) * (-intg(i  ,j  ,k-1,i_S_x)
                         -intg(i  ,j  ,k-1,i_S_y)
                         +intg(i  ,j  ,k-1,i_S_z))
         + Real(0.5) * (  intg(i  ,j  ,k-1,i_S_x_y)
                         -intg(i  ,j  ,k-1,i_S_x_z)
                         -intg(i  ,j  ,k-1,i_S_y_z))
            +          (  intg(i  ,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k-1) *
         ( Real(0.125) * vfrac(i-1,j  ,k-1)
         + Real(0.25) * ( intg(i-1,j  ,k-1,i_S_x)
                         -intg(i-1,j  ,k-1,i_S_y)
                         +intg(i-1,j  ,k-1,i_S_z))
         + Real(0.5) * ( -intg(i-1,j  ,k-1,i_S_x_y)
                         +intg(i-1,j  ,k-1,i_S_x_z)
                         -intg(i-1,j  ,k-1,i_S_y_z))
            +          ( -intg(i-1,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k-1) *
         ( Real(0.125) * vfrac(i  ,j-1,k-1)
         + Real(0.25) * (-intg(i  ,j-1,k-1,i_S_x)
                         +intg(i  ,j-1,k-1,i_S_y)
                         +intg(i  ,j-1,k-1,i_S_z))
         + Real(0.5) * ( -intg(i  ,j-1,k-1,i_S_x_y)
                         -intg(i  ,j-1,k-1,i_S_x_z)
                         +intg(i  ,j-1,k-1,i_S_y_z))
            +          ( -intg(i  ,j-1,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k-1) *
         ( Real(0.125) * vfrac(i-1,j-1,k-1)
         + Real(0.25) * ( intg(i-1,j-1,k-1,i_S_x)
                         +intg(i-1,j-1,k-1,i_S_y)
                         +intg(i-1,j-1,k-1,i_S_z))
         + Real(0.5) * (  intg(i-1,j-1,k-1,i_S_x_y)
                         +intg(i-1,j-1,k-1,i_S_x_z)
                         +intg(i-1,j-1,k-1,i_S_y_z))
             +         (  intg(i-1,j-1,k-1,i_S_xyz)));
    } else {
        return Real(0.);
    }
}

#endif

}
#endif
