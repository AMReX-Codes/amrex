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

template <typename T>
inline void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                             GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                             GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    if (domain.strictly_contains(vbx)) return;

    int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    Box const& sbox = amrex::grow(vbx,1);

    Box xlo_face = domain;
    xlo_face.setSmall(0,dlo.x-1);
    xlo_face.setBig  (0,dlo.x-1);
    xlo_face &= sbox;
    Box ylo_face = domain;
    ylo_face.setSmall(1,dlo.y-1);
    ylo_face.setBig  (1,dlo.y-1);
    ylo_face &= sbox;
    Box zlo_face = domain;
    zlo_face.setSmall(2,dlo.z-1);
    zlo_face.setBig  (2,dlo.z-1);
    zlo_face &= sbox;
    int xoffset = vbx.length(0)+1;
    int yoffset = vbx.length(1)+1;
    int zoffset = vbx.length(2)+1;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
    xlo_face, txbxlo,
    {
        auto lo = amrex::lbound(txbxlo);
        auto hi = amrex::ubound(txbxlo);
        if (lo.x == dlo.x-1 and bflo[0]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dlo.x-1,j,k) = a(dlo.x+offset,j,k);
            }}
        }
        if (lo.x+xoffset == dhi.x+1 and bfhi[0]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dhi.x+1,j,k) = a(dhi.x-offset,j,k);
            }}
        }
    },
    ylo_face, tybxlo,
    {
        auto lo = amrex::lbound(tybxlo);
        auto hi = amrex::ubound(tybxlo);
        if (lo.y == dlo.y-1 and bflo[1]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dlo.y-1,k) = a(i,dlo.y+offset,k);
            }}
        }
        if (lo.y+yoffset == dhi.y+1 and bfhi[1]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dhi.y+1,k) = a(i,dhi.y-offset,k);
            }}
        }
    },
    zlo_face, tzbxlo,
    {
        auto lo = amrex::lbound(tzbxlo);
        auto hi = amrex::ubound(tzbxlo);
        if (lo.z == dlo.z-1 and bflo[2]) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,j,dlo.z-1) = a(i,j,dlo.z+offset);
            }}
        }
        if (lo.z+zoffset == dhi.z+1 and bfhi[2]) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,j,dhi.z+1) = a(i,j,dhi.z-offset);
            }}
        }
    });

    const auto lo = amrex::lbound(sbox);
    const auto hi = amrex::ubound(sbox);

    AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
    {
        switch (iedge) {
        case 0: {
            // xlo & ylo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1) {
                if (bflo[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dlo.y-1,k) = a(dlo.x+offset,dlo.y-1,k);
                    }
                } else if (bflo[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dlo.y-1,k) = a(dlo.x-1,dlo.y+offset,k);
                    }
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1) {
                if (bfhi[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dlo.y-1,k) = a(dhi.x-offset,dlo.y-1,k);
                    }
                } else if (bflo[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dlo.y-1,k) = a(dhi.x+1,dlo.y+offset,k);
                    }
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1) {
                if (bflo[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dhi.y+1,k) = a(dlo.x+offset,dhi.y+1,k);
                    }
                } else if (bfhi[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dhi.y+1,k) = a(dlo.x-1,dhi.y-offset,k);
                    }
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1) {
                if (bfhi[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dhi.y+1,k) = a(dhi.x-offset,dhi.y+1,k);
                    }
                } else if (bfhi[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dhi.y+1,k) = a(dhi.x+1,dhi.y-offset,k);
                    }
                }
            }
            break;
        }
        case 4: {
            // xlo & zlo
            if (lo.x == dlo.x-1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dlo.z-1) = a(dlo.x+offset,j,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dlo.z-1) = a(dlo.x-1,j,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 5: {
            // xhi & zlo
            if (hi.x == dhi.x+1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dlo.z-1) = a(dhi.x-offset,j,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dlo.z-1) = a(dhi.x+1,j,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 6: {
            // xlo & zhi
            if (lo.x == dlo.x-1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dhi.z+1) = a(dlo.x+offset,j,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dhi.z+1) = a(dlo.x-1,j,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 7: {
            // xhi & zhi
            if (hi.x == dhi.x+1 and hi.z == dhi.z+1) {
                if (bfhi[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dhi.z+1) = a(dhi.x-offset,j,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dhi.z+1) = a(dhi.x+1,j,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 8: {
            // ylo & zlo
            if (lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bflo[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dlo.z-1) = a(i,dlo.y+offset,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dlo.z-1) = a(i,dlo.y-1,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 9: {
            // yhi & zlo
            if (hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bfhi[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dlo.z-1) = a(i,dhi.y-offset,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dlo.z-1) = a(i,dhi.y+1,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 10: {
            // ylo & zhi
            if (lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bflo[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dhi.z+1) = a(i,dlo.y+offset,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dhi.z+1) = a(i,dlo.y-1,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 11: {
            // yhi & zhi
            if (hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bfhi[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dhi.z+1) = a(i,dhi.y-offset,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dhi.z+1) = a(i,dhi.y+1,dhi.z-offset);
                    }
                }
            }
            break;
        }
        default: {}
        }
    });

    AMREX_HOST_DEVICE_FOR_1D ( 8, icorner,
    {
        switch (icorner) {
        case 0: {
            // xlo & ylo & zlo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x+offset,dlo.y-1,dlo.z-1);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x-1,dlo.y+offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x-1,dlo.y-1,dlo.z+offset);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo & zlo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x-offset,dlo.y-1,dlo.z-1);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x+1,dlo.y+offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x+1,dlo.y-1,dlo.z+offset);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi & zlo
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x+offset,dhi.y+1,dlo.z-1);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x-1,dhi.y-offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x-1,dhi.y+1,dlo.z+offset);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi & zlo
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x-offset,dhi.y+1,dlo.z-1);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x+1,dhi.y-offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x+1,dhi.y+1,dlo.z+offset);
                }
            }
            break;
        }
        case 4: {
            // xlo & ylo & zhi
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x+offset,dlo.y-1,dhi.z+1);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x-1,dlo.y+offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x-1,dlo.y-1,dhi.z-offset);
                }
            }
            break;
        }
        case 5: {
            // xhi & ylo & zhi
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x-offset,dlo.y-1,dhi.z+1);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x+1,dlo.y+offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x+1,dlo.y-1,dhi.z-offset);
                }
            }
            break;
        }
        case 6: {
            // xlo & yhi & zhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x+offset,dhi.y+1,dhi.z+1);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x-1,dhi.y-offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x-1,dhi.y+1,dhi.z-offset);
                }
            }
            break;
        }
        case 7: {
            // xhi & yhi & zhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bfhi[0] ) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x-offset,dhi.y+1,dhi.z+1);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x+1,dhi.y-offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x+1,dhi.y+1,dhi.z-offset);
                }
            }
            break;
        }
        default: {}
        }
    });
}

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
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc+1,kc))/(w1+w2);
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_z (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w2 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc,kc+1))/(w1+w2);
    };

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
    };

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
    };
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
    };

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
    };
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
        (  std::abs(sten(i-1,j,k,ist_p00)) + std::abs(sten(i,j,k,ist_p00))
         + std::abs(sten(i,j-1,k,ist_0p0)) + std::abs(sten(i,j,k,ist_0p0))
         + std::abs(sten(i,j,k-1,ist_00p)) + std::abs(sten(i,j,k,ist_00p))
         + std::abs(sten(i-1,j-1,k,ist_pp0)) + std::abs(sten(i,j-1,k,ist_pp0))
         + std::abs(sten(i-1,j,k,ist_pp0)) + std::abs(sten(i,j,k,ist_pp0))
         + std::abs(sten(i-1,j,k-1,ist_p0p)) + std::abs(sten(i,j,k-1,ist_p0p))
         + std::abs(sten(i-1,j,k,ist_p0p)) + std::abs(sten(i,j,k,ist_p0p))
         + std::abs(sten(i,j-1,k-1,ist_0pp)) + std::abs(sten(i,j,k-1,ist_0pp))
         + std::abs(sten(i,j-1,k,ist_0pp)) + std::abs(sten(i,j,k,ist_0pp))
         + std::abs(sten(i-1,j-1,k-1,ist_ppp)) + std::abs(sten(i,j-1,k-1,ist_ppp))
         + std::abs(sten(i-1,j,k-1,ist_ppp)) + std::abs(sten(i,j,k-1,ist_ppp))
         + std::abs(sten(i-1,j-1,k,ist_ppp)) + std::abs(sten(i,j-1,k,ist_ppp))
         + std::abs(sten(i-1,j,k,ist_ppp)) + std::abs(sten(i,j,k,ist_ppp)) + eps);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_sten (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
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
    });
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
            Real w1 = std::abs(sten(i,j,k-1,ist_00p));
            Real w2 = std::abs(sten(i,j,k  ,ist_00p));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc,kc+1));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc,kc+1)) / (w1+w2);
            }
        } else if (ieven and keven) {
            Real w1 = std::abs(sten(i,j-1,k,ist_0p0));
            Real w2 = std::abs(sten(i,j  ,k,ist_0p0));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc+1,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc+1,kc)) / (w1+w2);
            }
        } else if (jeven and keven) {
            Real w1 = std::abs(sten(i-1,j,k,ist_p00));
            Real w2 = std::abs(sten(i  ,j,k,ist_p00));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic+1,jc,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic+1,jc,kc)) / (w1+w2);
            }
        } else if (ieven) {
            Real w1m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j-1,k  ,ist_0pp)) + eps);
            Real w1p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i,j  ,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k-1,ist_0pp)) + eps);
            Real w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i,j-1,k  ,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real wmm = std::abs(sten(i,j-1,k-1,ist_0pp)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i,j  ,k-1,ist_0pp)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i,j-1,k  ,ist_0pp)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i,j  ,k  ,ist_0pp)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic,jc+1,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic,jc+1,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (jeven) {
            Real w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i-1,j,k  ,ist_p0p)) + eps);
            Real w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k-1,ist_p0p)) + eps);
            Real w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i-1,j,k  ,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real wmm = std::abs(sten(i-1,j,k-1,ist_p0p)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i  ,j,k-1,ist_p0p)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i-1,j,k  ,ist_p0p)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i  ,j,k  ,ist_p0p)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic+1,jc,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (keven) {
            Real w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i-1,j  ,k,ist_pp0)) + eps);
            Real w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real w2m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j-1,k,ist_pp0)) + eps);
            Real w2p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i-1,j  ,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real wmm = std::abs(sten(i-1,j-1,k,ist_pp0)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i  ,j-1,k,ist_pp0)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i-1,j  ,k,ist_pp0)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i  ,j  ,k,ist_pp0)) * (1.0 + w1p + w2p);
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

            Real wtmp = std::abs(sten(i-1,j,k,ist_p00)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i-1,j  ,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p00)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp))
                + std::abs(sten(i,j  ,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0p0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i  ,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0p0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp))
                + std::abs(sten(i  ,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i  ,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp))
                + std::abs(sten(i  ,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j-1,k,ist_pp0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmmp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_pp0)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_pp0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_pp0)) /
                ( std::abs(sten(i,j,k-1,ist_ppp))
                + std::abs(sten(i,j,k  ,ist_ppp)) + eps);
            wppm += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j,k-1,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_p0p)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p0p)) /
                ( std::abs(sten(i,j-1,k,ist_ppp))
                + std::abs(sten(i,j  ,k,ist_ppp)) + eps);
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0pp)) /
                ( std::abs(sten(i-1,j,k,ist_ppp))
                + std::abs(sten(i  ,j,k,ist_ppp)) + eps);
            wmpp += wtmp;
            wppp += wtmp;

            wmmm *= std::abs(sten(i-1,j-1,k-1,ist_ppp));
            wpmm *= std::abs(sten(i  ,j-1,k-1,ist_ppp));
            wmpm *= std::abs(sten(i-1,j  ,k-1,ist_ppp));
            wppm *= std::abs(sten(i  ,j  ,k-1,ist_ppp));
            wmmp *= std::abs(sten(i-1,j-1,k  ,ist_ppp));
            wpmp *= std::abs(sten(i  ,j-1,k  ,ist_ppp));
            wmpp *= std::abs(sten(i-1,j  ,k  ,ist_ppp));
            wppp *= std::abs(sten(i  ,j  ,k  ,ist_ppp));
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

        Real sten_lo = std::abs(sten(ii-2,jj,kk,ist_p00));
        Real sten_hi = std::abs(sten(ii-1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii-1,jj,kk);
        } else {
            cv += fine(ii-1,jj,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii+1,jj,kk)
        // ************************************

        sten_lo = std::abs(sten(ii  ,jj,kk,ist_p00));
        sten_hi = std::abs(sten(ii+1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii+1,jj,kk);
        } else {
            cv += fine(ii+1,jj,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj-1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj-2,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj-1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj-1,kk);
        } else {
            cv += fine(ii,jj-1,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj+1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj  ,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj+1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj+1,kk);
        } else {
            cv += fine(ii,jj+1,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk-1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk-2,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk-1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk-1);
        } else {
            cv += fine(ii,jj,kk-1)*sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk+1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk  ,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk+1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk+1);
        } else {
            cv += fine(ii,jj,kk+1)*sten_lo  / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii-1,jj-1,kk)
        // ************************************

        // keven
        Real w1m = std::abs(sten(ii-2,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-2,jj-1,kk,ist_pp0)) + eps);
        Real w1p = std::abs(sten(ii-1,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-1,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real w2m = std::abs(sten(ii-1,jj-2,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-2,kk,ist_pp0)) + eps);
        Real w2p = std::abs(sten(ii-1,jj-1,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-1,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real wmm = std::abs(sten(ii-2,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        Real wpm = std::abs(sten(ii-1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        Real wmp = std::abs(sten(ii-2,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        Real wpp = std::abs(sten(ii-1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj-1,kk)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj-1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj-1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii  ,jj-1,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii+1,jj-1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii+1,jj-2,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-2,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii+1,jj-1,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii  ,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj-1,kk)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii-2,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-2,jj+1,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii-1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii-1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj  ,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii-1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj+1,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii-2,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj+1,kk)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj+1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii+1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj  ,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii+1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii+1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj+1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii  ,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj+1,kk)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk-1)
        // ************************************

        // jeven
        w1m = std::abs(sten(ii-2,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-2,jj,kk-1,ist_p0p)) + eps);
        w1p = std::abs(sten(ii-1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        w2m = std::abs(sten(ii-1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-2,ist_p0p)) + eps);
        w2p = std::abs(sten(ii-1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii-2,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk-1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk-1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        w1p = std::abs(sten(ii+1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii+1,jj,kk-1,ist_p0p)) + eps);
        w2m = std::abs(sten(ii+1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-2,ist_p0p)) + eps);
        w2p = std::abs(sten(ii+1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii  ,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii-2,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-2,jj,kk  ,ist_p0p)) + eps);
        w1p = std::abs(sten(ii-1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2m = std::abs(sten(ii-1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii-2,jj,kk  ,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2p = std::abs(sten(ii-1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk+1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii-2,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk+1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w1p = std::abs(sten(ii+1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
              +std::abs(sten(ii+1,jj,kk  ,ist_p0p)) + eps);
        w2m = std::abs(sten(ii+1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii+1,jj,kk  ,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w2p = std::abs(sten(ii+1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk+1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii  ,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk-1)
        // ************************************

        // ieven
        w1m = std::abs(sten(ii,jj-2,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-2,kk-1,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj-1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-2,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj-1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj-1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj-2,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk-1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj+1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj+1,kk-1,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj+1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-2,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj+1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj  ,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj-2,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-2,kk  ,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj-1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj-1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj-2,kk  ,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj-1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk+1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj-2,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj+1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj+1,kk  ,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj+1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj+1,kk  ,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj+1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk+1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj  ,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine at corners
        // ************************************

        Real wmmm = 1.0
            + std::abs(sten(ii  ,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp)) + eps);
        wmmm *= std::abs(sten(ii,jj,kk,ist_ppp));
        cv += wmmm*fine(ii+1,jj+1,kk+1)*sten(ii+1,jj+1,kk+1,ist_inv);

        Real wpmm = 1.0
            + std::abs(sten(ii-1,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp)) + eps);
        wpmm *= std::abs(sten(ii-1,jj,kk,ist_ppp));
        cv += wpmm*fine(ii-1,jj+1,kk+1)*sten(ii-1,jj+1,kk+1,ist_inv);

        Real wmpm = 1.0
            + std::abs(sten(ii  ,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps);
        wmpm *= std::abs(sten(ii  ,jj-1,kk  ,ist_ppp));
        cv += wmpm*fine(ii+1,jj-1,kk+1)*sten(ii+1,jj-1,kk+1,ist_inv);

        Real wppm = 1.0
            + std::abs(sten(ii-1,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps);
        wppm *= std::abs(sten(ii-1,jj-1,kk  ,ist_ppp));
        cv += wppm*fine(ii-1,jj-1,kk+1)*sten(ii-1,jj-1,kk+1,ist_inv);

        Real wmmp = 1.0
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps);
        wmmp *= std::abs(sten(ii  ,jj  ,kk-1,ist_ppp));
        cv += wmmp*fine(ii+1,jj+1,kk-1)*sten(ii+1,jj+1,kk-1,ist_inv);

        Real wpmp = 1.0
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps);
        wpmp *= std::abs(sten(ii-1,jj  ,kk-1,ist_ppp));
        cv += wpmp*fine(ii-1,jj+1,kk-1)*sten(ii-1,jj+1,kk-1,ist_inv);

        Real wmpp = 1.0
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps);
        wmpp *= std::abs(sten(ii  ,jj-1,kk-1,ist_ppp));
        cv += wmpp*fine(ii+1,jj-1,kk-1)*sten(ii+1,jj-1,kk-1,ist_inv);

        Real wppp = 1.0
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps);
        wppp *= std::abs(sten(ii-1,jj-1,kk-1,ist_ppp));
        cv += wppp*fine(ii-1,jj-1,kk-1)*sten(ii-1,jj-1,kk-1,ist_inv);

        crse(i,j,k) = cv * 0.125;
    }
}

}
#endif
