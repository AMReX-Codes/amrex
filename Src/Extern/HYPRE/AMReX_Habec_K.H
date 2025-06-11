#ifndef AMREX_Habec_K_H_
#define AMREX_Habec_K_H_
#include <AMReX_Config.H>

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex::detail {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void comp_bf (Real& bf1, Real& bf2, Real sb, Real h, int bct, Real bcl, int bho)
    {
        Real fac = sb / (h*h);
        if (bct == AMREX_LO_DIRICHLET) {
            Real h2 = Real(0.5)*h;
            if (bho >= 1) {
                Real h3 = Real(3.0)*h2;
                bf1 = fac*((h3 - bcl)/(bcl + h2) - Real(1.0));
                bf2 = fac*(bcl - h2)/(bcl + h3);
            } else {
                bf1 = fac*( h/(bcl + h2) - Real(1.0));
                bf2 = Real(0.0);
            }
        } else if (bct == AMREX_LO_NEUMANN) {
            bf1 = -fac;
            bf2 = Real(0.0);
        } else {
            bf1 = bf2 = Real(0.0);
            AMREX_ASSERT_WITH_MESSAGE(false,"hpmat: unsupported boundary type");
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void comp_bflo (Real& bf1, Real& bf2, Real& bflo, Real sb, Real h, int bct, Real bcl, int bho)
    {
        Real fac = sb / (h*h);
        if (bct == AMREX_LO_DIRICHLET) {
            Real h2 = Real(0.5)*h;
            bflo = fac * ( h / (bcl + h2) - Real(1.0));
            if (bho >= 1) {
                Real h3 = Real(3.0)*h2;
                bf1 = fac*((h3 - bcl)/(bcl + h2) - Real(1.0));
                bf2 = fac*(bcl - h2)/(bcl + h3);
            } else {
                bf1 = bflo;
                bf2 = Real(0.0);
            }
        } else if (bct == AMREX_LO_NEUMANN) {
            bflo = -fac;
            bf1 = -fac;
            bf2 = Real(0.0);
        } else {
            bf1 = bf2 = bflo = Real(0.0);
            AMREX_ASSERT_WITH_MESSAGE(false,"hpmat: unsupported boundary type");
        }
    }
}

#if (AMREX_SPACEDIM == 2)
#include <AMReX_Habec_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_Habec_3D_K.H>
#endif

#endif
