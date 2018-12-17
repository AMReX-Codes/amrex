#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
inline
void
amrex_cellconslin_unlim_slopes (Box const& bx, FArrayBox& slopesfab,
                                FArrayBox const& ufab, const int icomp, const int ncomp,
                                BCRec const* AMREX_RESTRICT bcr)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto slopes = slopesfab.view(lo,0);
    const auto u = ufab.view(lo,icomp);
    const auto& sbox = slopesfab.box();
    const auto slo = amrex::lbound(sbox);
    const auto shi = amrex::ubound(sbox);
    const auto slen = amrex::length(sbox);

    for             (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n        ) = 0.5*(u(i+1,j,k,n)-u(i-1,j,k,n));
                    slopes(i,j,k,n+  ncomp) = 0.5*(u(i,j+1,k,n)-u(i,j-1,k,n));
                    slopes(i,j,k,n+2*ncomp) = 0.5*(u(i,j,k+1,n)-u(i,j,k-1,n));
                }
            }
        }

        const BCRec& bc = bcr[n];

        if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
        {
            const int i = 0;
            if (slen.x >= 2) {
                for         (int k = 0; k < len.z; ++k) {
                    for     (int j = 0; j < len.y; ++j) {
                        slopes(i,j,k,n) = -(16./15.)*u(i-1,j,k,n) + 0.5*u(i,j,k,n)
                            + (2./3.)*u(i+1,j,k,n) - 0.1*u(i+2,j,k,n);
                    }
                }
            } else {
                for         (int k = 0; k < len.z; ++k) {
                    for     (int j = 0; j < len.y; ++j) {
                        slopes(i,j,k,n) = 0.25*(u(i+1,j,k,n)+5.*u(i,j,k,n)-6.*u(i-1,j,k,n));
                    }
                }
            }
        }

        if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
        {
            const int i = len.x-1;
            if (slen.x >= 2) {
                for         (int k = 0; k < len.z; ++k) {
                    for     (int j = 0; j < len.y; ++j) {
                        slopes(i,j,k,n) = (16./15.)*u(i+1,j,k,n) - 0.5*u(i,j,k,n)
                            - (2./3.)*u(i-1,j,k,n) + 0.1*u(i-2,j,k,n);
                    }
                }                
            } else {
                for         (int k = 0; k < len.z; ++k) {
                    for     (int j = 0; j < len.y; ++j) {
                        slopes(i,j,k,n) = -0.25*(u(i-1,j,k,n)+5.*u(i,j,k,n)-6.*u(i+1,j,k,n));
                    }
                }
            }
        }

        if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
        {
            const int j = 0;
            if (slen.y >= 2) {
                for     (int k = 0; k < len.z; ++k) {
                    for (int i = 0; i < len.z; ++i) {
                        slopes(i,j,k,n+ncomp) = -(16./15.)*u(i,j-1,k,n) + 0.5*u(i,j,k,n)
                            + (2./3.)*u(i,j+1,k,n) - 0.1*u(i,j+2,k,n);
                    }
                }
            } else {
                for     (int k = 0; k < len.z; ++k) {
                    for (int i = 0; i < len.z; ++i) {
                        slopes(i,j,k,n+ncomp) = 0.25*(u(i,j+1,k,n)+5.*u(i,j,k,n)-6.*u(i,j-1,k,n));
                    }
                }
            }
        }

        if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
        {
            const int j = len.y-1;
            if (slen.y >= 2) {
                for     (int k = 0; k < len.z; ++k) {
                    for (int i = 0; i < len.z; ++i) {
                        slopes(i,j,k,n+ncomp) = (16./15.)*u(i,j+1,k,n) - 0.5*u(i,j,k,n)
                            - (2./3.)*u(i,j-1,k,n) + 0.1*u(i,j-2,k,n);
                    }
                }
            } else {
                for     (int k = 0; k < len.z; ++k) {
                    for (int i = 0; i < len.z; ++i) {
                        slopes(i,j,k,n+ncomp) = -0.25*(u(i,j-1,k,n)+5.*u(i,j,k,n)-6.*u(i,j+1,k,n));
                    }
                }
            }
        }

        if (lo.z == slo.z && (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::hoextrap))
        {
            const int k = 0;
            if (slen.z >= 2) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        slopes(i,j,k,n+2*ncomp) = -(16./15.)*u(i,j,k-1,n) + 0.5*u(i,j,k,n)
                            + (2./3.)*u(i,j,k+1,n) - 0.1*u(i,j,k+2,n);
                    }
                }
            } else {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        slopes(i,j,k,n+2*ncomp) = 0.25*(u(i,j,k+1,n)+5.*u(i,j,k,n)-6.*u(i,j,k-1,n));
                    }
                }
            }
        }

        if (hi.z == shi.z && (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::hoextrap))
        {
            const int k = len.z-1;
            if (slen.z >= 2) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        slopes(i,j,k,n+2*ncomp) = (16./15.)*u(i,j,k+1,n) - 0.5*u(i,j,k,n)
                            - (2./3.)*u(i,j,k-1,n) + 0.1*u(i,j,k-2,n);
                    }
                }
            } else {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        slopes(i,j,k,n+2*ncomp) = -0.25*(u(i,j,k-1,n)+5.*u(i,j,k,n)-6.*u(i,j,k+1,n));
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void
amrex_cellconslin_unlim_interp (Box const& bx, FArrayBox& finefab, const int fcomp, const int ncomp,
                                FArrayBox const& slopesfab, FArrayBox const& crsefab, const int ccomp,
                                IntVect const& ratio)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto fine = finefab.view(lo,fcomp);

    const auto clo = amrex::coarsen(lo,ratio);
    const auto slopes = slopesfab.view(clo);
    const auto crse = crsefab.view(clo,ccomp);

    const Real xrinv = 1./ratio[0];
    const Real yrinv = 1./ratio[1];
    const Real zrinv = 1./ratio[2];

    const Real xtmp = 0.5*(xrinv-1.0);
    const Real ytmp = 0.5*(yrinv-1.0);
    const Real ztmp = 0.5*(xrinv-1.0);


    for (int k = 0; k < len.z; ++k) {
        const int kk = k + lo.z;
        const int kkc = amrex::coarsen(kk, ratio[2]);
        const int koff = kk - kkc*ratio[2];
        const Real zoff = koff*zrinv + ztmp;
        const int kc = kkc - clo.z;
        for (int j = 0; j < len.y; ++j) {
            const int jj = j + lo.y;
            const int jjc = amrex::coarsen(jj, ratio[1]);
            const int joff = jj - jjc*ratio[1];
            const Real yoff = joff*yrinv + ytmp;
            const int jc = jjc - clo.y;
#ifndef __CUDA_ARCH__
            for (int n = 0; n < ncomp; ++n) {
#endif
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const int ii = i + lo.x;
                const int iic = amrex::coarsen(ii, ratio[0]);
                const int ioff = ii - iic*ratio[0];
                const Real xoff = ioff*xrinv + xtmp;
                const int ic = iic - clo.x;
#ifdef __CUDA_ARCH__
                for (int n = 0; n < ncomp; ++n) {
#endif
                    fine(i,j,k,n) = crse(ic,jc,kc,n)
                        + xoff * slopes(ic,jc,kc,n)
                        + yoff * slopes(ic,jc,kc,n+ncomp)
                        + zoff * slopes(ic,jc,kc,n+2*ncomp);
                }
            }
        }
    }
}

}

#endif
