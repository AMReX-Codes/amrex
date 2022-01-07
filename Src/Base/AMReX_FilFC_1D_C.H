#ifndef AMREX_FILFC_1D_C_H_
#define AMREX_FILFC_1D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

struct FilfcFace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const IntVect& iv, Array4<Real> const& q,
                     const int dcomp, const int numcomp,
                     Box const& domain_box, const BCRec* bcr,
                     const int bcomp) const noexcept
    {
        const int i = iv[0];

        // Domain box is indexed according to the face currently treating
        const auto& domain_lo = domain_box.loVect();
        const auto& domain_hi = domain_box.hiVect();
        const int ilo = domain_lo[0];
        const int ihi = domain_hi[0];

        for (int n = dcomp; n < numcomp+dcomp; ++n)
        {
            const BCRec& bc = bcr[bcomp+n-dcomp];

            if (i == ilo)
            {
                // Enforce reflect_odd on the domain face
                if (bc.lo(0) == BCType::reflect_odd) {
                    q(i,0,0,n) = 0.0;
                }
            }
            else if (i < ilo)
            {
                switch (bc.lo(0)) {
                case (BCType::foextrap):
                {
                    q(i,0,0,n) = q(ilo,0,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i < ilo - 1)
                    {
                        q(i,0,0,n) = q(ilo,0,0,n);
                    }
                    // i == ilo-1
                    else
                    {
                        q(i,0,0,n) = Real(2.0)*q(i+1,0,0,n) - q(i+2,0,0,n);
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,0,0,n) = q(2*ilo-i,0,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,0,0,n) = -q(2*ilo-i,0,0,n);
                    break;
                }
                }
            }
            else if (i == ihi)
            {
                // Enforce reflect_odd on the domain face
                if (bc.hi(0) == BCType::reflect_odd) {
                    q(i,0,0,n) = 0.0;
                }
            }
            else if (i > ihi)
            {
                switch (bc.hi(0)) {
                case (BCType::foextrap):
                {
                    q(i,0,0,n) = q(ihi,0,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i > ihi + 1)
                    {
                        q(i,0,0,n) = q(ihi,0,0,n);
                    }
                    // i == ihi+1
                    else
                    {
                        q(i,0,0,n) = Real(2.0)*q(i-1,0,0,n) - q(i-2,0,0,n);
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,0,0,n) = q(2*ihi-i,0,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,0,0,n) = -q(2*ihi-i,0,0,n);
                    break;
                }
                }
            }
        }
    }
};

}

#endif
