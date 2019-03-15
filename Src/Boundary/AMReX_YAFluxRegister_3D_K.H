#ifndef AMREX_YAFLUXREGISTER_3D_K_H_
#define AMREX_YAFLUXREGISTER_3D_K_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void yafluxreg_crseadd (Box const& bx, Array4<Real> const& d, Array4<int const> const& flag,
                        FArrayBox const& fxfab, FArrayBox const& fyfab, FArrayBox const& fzfab,
                        Real dtdx, Real dtdy, Real dtdz, int nc) noexcept
{
    auto const fx = fxfab.array();
    auto const fy = fyfab.array();
    auto const fz = fzfab.array();

    auto const lo = amrex::lbound(bx);
    auto const hi = amrex::ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (flag(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell)
        {
            if (flag(i-1,j,k) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) -= dtdx*fx(i,j,k,n);
                }
            } else if (flag(i+1,j,k) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) += dtdx*fx(i+1,j,k,n);
                }
            }

            if (flag(i,j-1,k) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) -= dtdy*fy(i,j,k,n);
                }
            } else if (flag(i,j+1,k) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) += dtdy*fy(i,j+1,k,n);
                }
            }

            if (flag(i,j,k-1) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) -= dtdz*fz(i,j,k,n);
                }
            } else if (flag(i,j,k+1) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,k,n) += dtdz*fz(i,j,k+1,n);
                }
            }
        }
    }}}
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void yafluxreg_fineadd (Box const& bx, Array4<Real> const& d, FArrayBox const& ffab,
                        Real dtdx, int nc, int dirside, Dim3 const& rr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto f = ffab.array();

    switch (dirside) {
    case 0 :
    {
        for (int n = 0; n < nc; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    const int i = lo.x;
                    const int ii = (i+1)*rr.x;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (    int koff = 0; koff < rr.z; ++koff) {
                        const     int kk =      k*rr.z  + koff;
                        for (int joff = 0; joff < rr.y; ++joff) {
                            const int jj =      j*rr.y  + joff;
                            Real tmp = -dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                            *dp += tmp;
                        }
                    }
                }
            }
        }
        break;
    }
    case 1 :
    {
        for (int n = 0; n < nc; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    const int i = lo.x;
                    const int ii = i*rr.x;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (int     koff = 0; koff < rr.z; ++koff) {
                        const     int kk =      k*rr.z  + koff;
                        for (int joff = 0; joff < rr.y; ++joff) {
                            const int jj =      j*rr.y  + joff;
                            Real tmp = dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                            *dp += tmp;
                        }
                    }
                }
            }
        }
        break;
    }
    case 2 :
    {
        for (int n = 0; n < nc; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int j = lo.y;
                    const int jj = (j+1)*rr.y;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (    int koff = 0; koff < rr.z; ++koff) {
                        const     int kk =      k*rr.z  + koff;
                        for (int ioff = 0; ioff < rr.x; ++ioff) {
                            const int ii =      i*rr.x  + ioff;
                            Real tmp = -dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                            *dp += tmp;
                        }
                    }
                }
            }
        }
        break;
    }
    case 3 :
    {
        for (int n = 0; n < nc; ++n) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int j = lo.y;
                    const int jj = j*rr.y;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (    int koff = 0; koff < rr.z; ++koff) {
                        const     int kk =      k*rr.z  + koff;
                        for (int ioff = 0; ioff < rr.x; ++ioff) {
                            const int ii =      i*rr.x  + ioff;
                            Real tmp = dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                            *dp += tmp;
                        }
                    }
                }
            }
        }
        break;
    }
    case 4 :
    {
        for (int n = 0; n < nc; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int k = lo.z;
                    const int kk = (k+1)*rr.z;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (    int joff = 0; joff < rr.y; ++joff) {
                        const     int jj =      j*rr.y  + joff;
                        for (int ioff = 0; ioff < rr.x; ++ioff) {
                            const int ii =      i*rr.x  + ioff;
                            Real tmp = -dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                            *dp += tmp;
                        }
                    }
                }
            }
        }
        break;
    }
    default:
    {
        for (int n = 0; n < nc; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int k = lo.z;
                    const int kk = k*rr.z;
                    Real* AMREX_RESTRICT dp = &(d(i,j,k,n));
                    for (    int joff = 0; joff < rr.y; ++joff) {
                        const     int jj =      j*rr.y  + joff;
                        for (int ioff = 0; ioff < rr.x; ++ioff) {
                            const int ii =      i*rr.x  + ioff;
                            Real tmp = dtdx*f(ii,jj,kk,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                                *dp += tmp;
                        }
                    }
                }
            }
        }
    }
    }
}

}
#endif
