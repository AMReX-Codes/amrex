#ifndef AMREX_FILCC_1D_C_H_
#define AMREX_FILCC_1D_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

AMREX_GPU_DEVICE
inline
void
filcc_cell (const IntVect& iv, Array4<Real> const& q,
            const int dcomp, const int numcomp,
            GeometryData const& geom, const Real /*time*/,
            const BCRec* bcr, const int bcomp,
            const int /*orig_comp*/) noexcept
{
    const int i = iv[0];

    const Box& domain_box = geom.Domain();
    const auto& domain_lo = domain_box.loVect();
    const auto& domain_hi = domain_box.hiVect();
    const int ilo = domain_lo[0];
    const int ihi = domain_hi[0];
    const int is = amrex::max(q.begin.x,ilo);
    const int ie = amrex::min(q.end.x-1,ihi);

    for (int n = dcomp; n < numcomp+dcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n-dcomp];

        if (i < ilo)
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
                else if (ilo+2 <= ie)
                {
                    q(i,0,0,n) = 0.125*(15.*q(i+1,0,0,n) - 10.*q(i+2,0,0,n) + 3.*q(i+3,0,0,n));
                }
                else
                {
                    q(i,0,0,n) = 0.5*(3.*q(i+1,0,0,n) - q(i+2,0,0,n));
                }
                break;
            }
            case (BCType::hoextrapcc):
            {
                q(i,0,0,n) = 2.*q(i+1,0,0,n) - q(i+2,0,0,n);
                break;
            }
            case (BCType::reflect_even):
            {
                q(i,0,0,n) = q(2*ilo-i-1,0,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(i,0,0,n) = -q(2*ilo-i-1,0,0,n);
                break;
            }
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
                else if (ihi-2 >= is)
                {
                    q(i,0,0,n) = 0.125*(15.*q(i-1,0,0,n) - 10.*q(i-2,0,0,n) + 3.*q(i-3,0,0,n));
                }
                else
                {
                    q(i,0,0,n) = 0.5*(3.*q(i-1,0,0,n) - q(i-2,0,0,n));
                }
                break;
            }
            case (BCType::hoextrapcc):
            {
                q(i,0,0,n) = 2.*q(i-1,0,0,n) - q(i-2,0,0,n);
                break;
            }
            case (BCType::reflect_even):
            {
                q(i,0,0,n) = q(2*ihi-i+1,0,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(i,0,0,n) = -q(2*ihi-i+1,0,0,n);
                break;
            }
            }
        }
    }
}

}

#endif

