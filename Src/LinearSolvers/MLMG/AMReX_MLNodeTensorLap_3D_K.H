#ifndef AMREX_MLNODETENSORLAP_3D_K_H_
#define AMREX_MLNODETENSORLAP_3D_K_H_
#include <AMReX_Config.H>

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
        if (i_is_odd && j_is_odd && k_is_odd) {
            // Fine node at center of cell
            fine(i,j,k) += (ts_interp_face_yz(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_yz(crse,ic+1,jc  ,kc  ) +
                            ts_interp_face_xz(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_xz(crse,ic  ,jc+1,kc  ) +
                            ts_interp_face_xy(crse,ic  ,jc  ,kc  ) +
                            ts_interp_face_xy(crse,ic  ,jc  ,kc+1)) * Real(1./6.);
        } else if (j_is_odd && k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += ts_interp_face_yz(crse,ic,jc,kc);
        } else if (i_is_odd && k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += ts_interp_face_xz(crse,ic,jc,kc);
        } else if (i_is_odd && j_is_odd) {
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
void mlndtslap_semi_interpadd (int i, int j, int k, Array4<Real> const& fine,
                               Array4<Real const> const& crse, Array4<int const> const& msk,
                               int semi_dir) noexcept
{
    if (!msk(i,j,k)) {
        if (semi_dir == 0) {
            int jc = amrex::coarsen(j,2);
            int kc = amrex::coarsen(k,2);
            bool j_is_odd = (jc*2 != j);
            bool k_is_odd = (kc*2 != k);
            if (j_is_odd && k_is_odd) {
                // Node on a Y-Z face
                fine(i,j,k) += ts_interp_face_yz(crse,i,jc,kc);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += ts_interp_line_y(crse,i,jc,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += ts_interp_line_z(crse,i,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(i,jc,kc);
            }
        } else if (semi_dir == 1) {
            int ic = amrex::coarsen(i,2);
            int kc = amrex::coarsen(k,2);
            bool i_is_odd = (ic*2 != i);
            bool k_is_odd = (kc*2 != k);
            if (i_is_odd && k_is_odd) {
                // Node on a Z-X face
                fine(i,j,k) += ts_interp_face_xz(crse,ic,j,kc);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += ts_interp_line_x(crse,ic,j,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += ts_interp_line_z(crse,ic,j,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,j,kc);
            }
        } else {
            int ic = amrex::coarsen(i,2);
            int jc = amrex::coarsen(j,2);
            bool i_is_odd = (ic*2 != i);
            bool j_is_odd = (jc*2 != j);
            if (i_is_odd && j_is_odd) {
                // Node on a X-Y face
                fine(i,j,k) += ts_interp_face_xy(crse,ic,jc,k);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += ts_interp_line_x(crse,ic,jc,k);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += ts_interp_line_y(crse,ic,jc,k);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,k);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                      Array4<int const> const& msk, GpuArray<Real,6> const& s) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        y(i,j,k) = s[0] * (x(i-1,j  ,k  ) + x(i+1,j  ,k  ))
            +      s[3] * (x(i  ,j-1,k  ) + x(i  ,j+1,k  ))
            +      s[5] * (x(i  ,j  ,k-1) + x(i  ,j  ,k+1))
            - Real(2.)*(s[0]+s[3]+s[5]) * x(i,j,k)
            + Real(0.5)*s[1] * (x(i-1,j-1,k  ) + x(i+1,j+1,k  ) - x(i-1,j+1,k  ) - x(i+1,j-1,k  ))
            + Real(0.5)*s[2] * (x(i-1,j  ,k-1) + x(i+1,j  ,k+1) - x(i-1,j  ,k+1) - x(i+1,j  ,k-1))
            + Real(0.5)*s[4] * (x(i  ,j-1,k-1) + x(i  ,j+1,k+1) - x(i  ,j-1,k+1) - x(i  ,j+1,k-1));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_gauss_seidel (int i, int j, int k, Array4<Real> const& sol,
                             Array4<Real const> const& rhs, Array4<int const> const& msk,
                             GpuArray<Real,6> const& s) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = 0.0;
    } else {
        constexpr Real omega = Real(1.25);
        Real s0 = Real(-2.)*(s[0]+s[3]+s[5]);
        Real Ax = s[0] * (sol(i-1,j  ,k  ) + sol(i+1,j  ,k  ))
            +     s[3] * (sol(i  ,j-1,k  ) + sol(i  ,j+1,k  ))
            +     s[5] * (sol(i  ,j  ,k-1) + sol(i  ,j  ,k+1))
            + s0 * sol(i,j,k)
            + Real(0.5)*s[1] * (sol(i-1,j-1,k  ) + sol(i+1,j+1,k  ) - sol(i-1,j+1,k  ) - sol(i+1,j-1,k  ))
            + Real(0.5)*s[2] * (sol(i-1,j  ,k-1) + sol(i+1,j  ,k+1) - sol(i-1,j  ,k+1) - sol(i+1,j  ,k-1))
            + Real(0.5)*s[4] * (sol(i  ,j-1,k-1) + sol(i  ,j+1,k+1) - sol(i  ,j-1,k+1) - sol(i  ,j+1,k-1));
        sol(i,j,k) += (rhs(i,j,k) - Ax) * (omega/s0);
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)

template <typename HypreInt, typename AtomicInt>
void mlndtslap_fill_ijmatrix_cpu (Box const& ndbx,
                                  Array4<AtomicInt const> const& gid,
                                  Array4<int const> const& lid,
                                  HypreInt* const ncols, HypreInt* const cols,
                                  Real* const mat, // NOLINT(readability-non-const-parameter)
                                  GpuArray<Real,6> const& s) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AtomicInt>::max();
    HypreInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            HypreInt nelems_old = nelems;

            cols[nelems] = gid(i,j,k);
            mat[nelems] = Real(-2.)*(s[0]+s[3]+s[5]);
            ++nelems;

            if (               gid(i  ,j-1,k-1) < gidmax) {
                cols[nelems] = gid(i  ,j-1,k-1);
                mat[nelems] = Real(0.5)*s[4];
                ++nelems;
            }

            if (               gid(i-1,j  ,k-1) < gidmax) {
                cols[nelems] = gid(i-1,j  ,k-1);
                mat[nelems] = Real(0.5)*s[2];
                ++nelems;
            }

            if (               gid(i  ,j  ,k-1) < gidmax) {
                cols[nelems] = gid(i  ,j  ,k-1);
                mat[nelems] = s[5];
                ++nelems;
            }

            if (               gid(i+1,j  ,k-1) < gidmax) {
                cols[nelems] = gid(i+1,j  ,k-1);
                mat[nelems] = Real(-0.5)*s[2];
                ++nelems;
            }

            if (               gid(i  ,j+1,k-1) < gidmax) {
                cols[nelems] = gid(i  ,j+1,k-1);
                mat[nelems] = Real(-0.5)*s[4];
                ++nelems;
            }

            if (               gid(i-1,j-1,k  ) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k  );
                mat[nelems] = Real(0.5)*s[1];
                ++nelems;
            }

            if (               gid(i  ,j-1,k  ) < gidmax) {
                cols[nelems] = gid(i  ,j-1,k  );
                mat[nelems] = s[3];
                ++nelems;
            }

            if (               gid(i+1,j-1,k  ) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k  );
                mat[nelems] = Real(-0.5)*s[1];
                ++nelems;
            }

            if (               gid(i-1,j  ,k  ) < gidmax) {
                cols[nelems] = gid(i-1,j  ,k  );
                mat[nelems] = s[0];
                ++nelems;
            }

            if (               gid(i+1,j  ,k  ) < gidmax) {
                cols[nelems] = gid(i+1,j  ,k  );
                mat[nelems] = s[0];
                ++nelems;
            }

            if (               gid(i-1,j+1,k  ) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k  );
                mat[nelems] = Real(-0.5)*s[1];
                ++nelems;
            }

            if (               gid(i  ,j+1,k  ) < gidmax) {
                cols[nelems] = gid(i  ,j+1,k  );
                mat[nelems] = s[3];
                ++nelems;
            }

            if (               gid(i+1,j+1,k  ) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k  );
                mat[nelems] = Real(0.5)*s[1];
                ++nelems;
            }

            if (               gid(i  ,j-1,k+1) < gidmax) {
                cols[nelems] = gid(i  ,j-1,k+1);
                mat[nelems] = Real(-0.5)*s[4];
                ++nelems;
            }

            if (               gid(i-1,j  ,k+1) < gidmax) {
                cols[nelems] = gid(i-1,j  ,k+1);
                mat[nelems] = Real(-0.5)*s[2];
                ++nelems;
            }

            if (               gid(i  ,j  ,k+1) < gidmax) {
                cols[nelems] = gid(i  ,j  ,k+1);
                mat[nelems] = s[5];
                ++nelems;
            }

            if (               gid(i+1,j  ,k+1) < gidmax) {
                cols[nelems] = gid(i+1,j  ,k+1);
                mat[nelems] = Real(0.5)*s[2];
                ++nelems;
            }

            if (               gid(i  ,j+1,k+1) < gidmax) {
                cols[nelems] = gid(i  ,j+1,k+1);
                mat[nelems] = Real(0.5)*s[4];
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
                                  HypreInt* const ncols, HypreInt* const cols,
                                  Real* const mat, // NOLINT(readability-non-const-parameter)
                                  GpuArray<Real,6> const& s) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        constexpr auto gidmax = std::numeric_limits<AtomicInt>::max();

        if (offset == 0) {
            cols[ps] = gid(i,j,k);
            mat[ps] = Real(-2.)*(s[0]+s[3]+s[5]);
            int nc = 1;
            if (gid(i-1,j-1,k-1) < gidmax) { ++nc; }
            if (gid(i  ,j-1,k-1) < gidmax) { ++nc; }
            if (gid(i+1,j-1,k-1) < gidmax) { ++nc; }
            if (gid(i-1,j  ,k-1) < gidmax) { ++nc; }
            if (gid(i  ,j  ,k-1) < gidmax) { ++nc; }
            if (gid(i+1,j  ,k-1) < gidmax) { ++nc; }
            if (gid(i-1,j+1,k-1) < gidmax) { ++nc; }
            if (gid(i  ,j+1,k-1) < gidmax) { ++nc; }
            if (gid(i+1,j+1,k-1) < gidmax) { ++nc; }
            if (gid(i-1,j-1,k  ) < gidmax) { ++nc; }
            if (gid(i  ,j-1,k  ) < gidmax) { ++nc; }
            if (gid(i+1,j-1,k  ) < gidmax) { ++nc; }
            if (gid(i-1,j  ,k  ) < gidmax) { ++nc; }
            if (gid(i+1,j  ,k  ) < gidmax) { ++nc; }
            if (gid(i-1,j+1,k  ) < gidmax) { ++nc; }
            if (gid(i  ,j+1,k  ) < gidmax) { ++nc; }
            if (gid(i+1,j+1,k  ) < gidmax) { ++nc; }
            if (gid(i-1,j-1,k+1) < gidmax) { ++nc; }
            if (gid(i  ,j-1,k+1) < gidmax) { ++nc; }
            if (gid(i+1,j-1,k+1) < gidmax) { ++nc; }
            if (gid(i-1,j  ,k+1) < gidmax) { ++nc; }
            if (gid(i  ,j  ,k+1) < gidmax) { ++nc; }
            if (gid(i+1,j  ,k+1) < gidmax) { ++nc; }
            if (gid(i-1,j+1,k+1) < gidmax) { ++nc; }
            if (gid(i  ,j+1,k+1) < gidmax) { ++nc; }
            if (gid(i+1,j+1,k+1) < gidmax) { ++nc; }
            ncols[lid(i,j,k)] = nc;
        }
        else if (offset == 1 && gid(i-1,j-1,k-1) < gidmax) {
            cols[ps] =          gid(i-1,j-1,k-1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 2 && gid(i  ,j-1,k-1) < gidmax) {
            cols[ps] =          gid(i  ,j-1,k-1);
            mat[ps] = Real(0.5)*s[4];
        }
        else if (offset == 3 && gid(i+1,j-1,k-1) < gidmax) {
            cols[ps] =          gid(i+1,j-1,k-1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 4 && gid(i-1,j  ,k-1) < gidmax) {
            cols[ps] =          gid(i-1,j  ,k-1);
            mat[ps] = Real(0.5)*s[2];
        }
        else if (offset == 5 && gid(i  ,j  ,k-1) < gidmax) {
            cols[ps] =          gid(i  ,j  ,k-1);
            mat[ps] = s[5];
        }
        else if (offset == 6 && gid(i+1,j  ,k-1) < gidmax) {
            cols[ps] =          gid(i+1,j  ,k-1);
            mat[ps] = Real(-0.5)*s[2];
        }
        else if (offset == 7 && gid(i-1,j+1,k-1) < gidmax) {
            cols[ps] =          gid(i-1,j+1,k-1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 8 && gid(i  ,j+1,k-1) < gidmax) {
            cols[ps] =          gid(i  ,j+1,k-1);
            mat[ps] = Real(-0.5)*s[4];
        }
        else if (offset == 9 && gid(i+1,j+1,k-1) < gidmax) {
            cols[ps] =          gid(i+1,j+1,k-1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 10 && gid(i-1,j-1,k  ) < gidmax) {
            cols[ps] =           gid(i-1,j-1,k  );
            mat[ps] = Real(0.5)*s[1];
        }
        else if (offset == 11 && gid(i  ,j-1,k  ) < gidmax) {
            cols[ps] =           gid(i  ,j-1,k  );
            mat[ps] = s[3];
        }
        else if (offset == 12 && gid(i+1,j-1,k  ) < gidmax) {
            cols[ps] =           gid(i+1,j-1,k  );
            mat[ps] = Real(-0.5)*s[1];
        }
        else if (offset == 13 && gid(i-1,j  ,k  ) < gidmax) {
            cols[ps] =           gid(i-1,j  ,k  );
            mat[ps] = s[0];
        }
        else if (offset == 14 && gid(i+1,j  ,k  ) < gidmax) {
            cols[ps] =           gid(i+1,j  ,k  );
            mat[ps] = s[0];
        }
        else if (offset == 15 && gid(i-1,j+1,k  ) < gidmax) {
            cols[ps] =           gid(i-1,j+1,k  );
            mat[ps] = Real(-0.5)*s[1];
        }
        else if (offset == 16 && gid(i  ,j+1,k  ) < gidmax) {
            cols[ps] =           gid(i  ,j+1,k  );
            mat[ps] = s[3];
        }
        else if (offset == 17 && gid(i+1,j+1,k  ) < gidmax) {
            cols[ps] =           gid(i+1,j+1,k  );
            mat[ps] = Real(0.5)*s[1];
        }
        else if (offset == 18 && gid(i-1,j-1,k+1) < gidmax) {
            cols[ps] =           gid(i-1,j-1,k+1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 19 && gid(i  ,j-1,k+1) < gidmax) {
            cols[ps] =           gid(i  ,j-1,k+1);
            mat[ps] = Real(-0.5)*s[4];
        }
        else if (offset == 20 && gid(i+1,j-1,k+1) < gidmax) {
            cols[ps] =           gid(i+1,j-1,k+1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 21 && gid(i-1,j  ,k+1) < gidmax) {
            cols[ps] =           gid(i-1,j  ,k+1);
            mat[ps] = Real(-0.5)*s[2];
        }
        else if (offset == 22 && gid(i  ,j  ,k+1) < gidmax) {
            cols[ps] =           gid(i  ,j  ,k+1);
            mat[ps] = s[5];
        }
        else if (offset == 23 && gid(i+1,j  ,k+1) < gidmax) {
            cols[ps] =           gid(i+1,j  ,k+1);
            mat[ps] = Real(0.5)*s[2];
        }
        else if (offset == 24 && gid(i-1,j+1,k+1) < gidmax) {
            cols[ps] =           gid(i-1,j+1,k+1);
            mat[ps] = Real(0.0);
        }
        else if (offset == 25 && gid(i  ,j+1,k+1) < gidmax) {
            cols[ps] =           gid(i  ,j+1,k+1);
            mat[ps] = Real(0.5)*s[4];
        }
        else if (offset == 26 && gid(i+1,j+1,k+1) < gidmax) {
            cols[ps] =           gid(i+1,j+1,k+1);
            mat[ps] = Real(0.0);
        }
    }
}

#endif

#endif

}

#endif
