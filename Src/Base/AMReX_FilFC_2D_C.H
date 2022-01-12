#ifndef AMREX_FILFC_2D_C_H_
#define AMREX_FILFC_2D_C_H_
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
        const int j = iv[1];

        // Domain box is indexed according to the face currently treating
        const IndexType &idxType = domain_box.ixType();
        const auto& domain_lo = domain_box.loVect();
        const auto& domain_hi = domain_box.hiVect();
        const int ilo = domain_lo[0];
        const int jlo = domain_lo[1];
        const int ihi = domain_hi[0];
        const int jhi = domain_hi[1];

        for (int n = dcomp; n < numcomp+dcomp; ++n)
        {
            const BCRec& bc = bcr[bcomp+n-dcomp];

            if (i == ilo)
            {
                // Enforce reflect_odd on the x domain face
                if ((bc.lo(0) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(0)) ) {
                    q(i,j,0,n) = 0.0;
                }
            }
            else if (i < ilo)
            {
                switch (bc.lo(0)) {
                case (BCType::foextrap):
                {
                    q(i,j,0,n) = q(ilo,j,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i < ilo - 1)
                    {
                        q(i,j,0,n) = q(ilo,j,0,n);
                    }
                    // i == ilo-1
                    else
                    {
                        q(i,j,0,n) = (idxType.nodeCentered(0)) ? Real(2.0)*q(i+1,j,0,n) - q(i+2,j,0,n)
                                                               : Real(0.5)*(Real(3.)*q(i+1,j,0,n) - q(i+2,j,0,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(0)) ? q(2*ilo-i,j,0,n)
                                                           : q(2*ilo-i-1,j,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(0)) ? -q(2*ilo-i,j,0,n)
                                                           : -q(2*ilo-i-1,j,0,n);
                    break;
                }
                }
            }
            else if (i == ihi)
            {
                // Enforce reflect_odd on the x domain face
                if ((bc.hi(0) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(0)) ) {
                    q(i,j,0,n) = 0.0;
                }
            }
            else if (i > ihi)
            {
                switch (bc.hi(0)) {
                case (BCType::foextrap):
                {
                    q(i,j,0,n) = q(ihi,j,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i > ihi + 1)
                    {
                        q(i,j,0,n) = q(ihi,j,0,n);
                    }
                    // i == ihi+1
                    else
                    {
                        q(i,j,0,n) = (idxType.nodeCentered(0)) ? Real(2.0)*q(i-1,j,0,n) - q(i-2,j,0,n)
                                                               : Real(0.5)*(Real(3.)*q(i-1,j,0,n) - q(i-2,j,0,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(0)) ? q(2*ihi-i,j,0,n)
                                                           : q(2*ihi-i+1,j,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(0)) ? -q(2*ihi-i,j,0,n)
                                                           : -q(2*ihi-i+1,j,0,n);
                    break;
                }
                }
            }

            if (j == jlo)
            {
                // Enforce reflect_odd on the y domain face
                if ((bc.lo(1) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(1)) ) {
                    q(i,j,0,n) = 0.0;
                }
            }
            else if (j < jlo)
            {
                switch (bc.lo(1)) {
                case (BCType::foextrap):
                {
                    q(i,j,0,n) = q(i,jlo,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (j < jlo - 1)
                    {
                        q(i,j,0,n) = q(i,jlo,0,n);
                    }
                    // j == jlo-1
                    else
                    {
                        q(i,j,0,n) = (idxType.nodeCentered(1)) ? Real(2.0)*q(i,j+1,0,n) - q(i,j+2,0,n)
                                                               : Real(0.5)*(Real(3.)*q(i,j+1,0,n) - q(i,j+2,0,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(1)) ? q(i,2*jlo-j,0,n)
                                                           : q(i,2*jlo-j-1,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(1)) ? -q(i,2*jlo-j,0,n)
                                                           : -q(i,2*jlo-j-1,0,n);
                    break;
                }
                }
            }
            else if (j == jhi)
            {
                // Enforce reflect_odd on the y domain face
                if ((bc.hi(1) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(1)) ) {
                    q(i,j,0,n) = 0.0;
                }
            }
            else if (j > jhi)
            {
                switch (bc.hi(1)) {
                case (BCType::foextrap):
                {
                    q(i,j,0,n) = q(i,jhi,0,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (j > jhi + 1)
                    {
                        q(i,j,0,n) = q(i,jhi,0,n);
                    }
                    // j == jhi+1
                    else
                    {
                        q(i,j,0,n) = (idxType.nodeCentered(1)) ? Real(2.0)*q(i,j-1,0,n) - q(i,j-2,0,n)
                                                               : Real(0.5)*(Real(3.)*q(i,j-1,0,n) - q(i,j-2,0,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(1)) ? q(i,2*jhi-j,0,n)
                                                           : q(i,2*jhi-j+1,0,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,0,n) = (idxType.nodeCentered(1)) ? -q(i,2*jhi-j,0,n)
                                                           : -q(i,2*jhi-j+1,0,n);
                    break;
                }
                }
            }
        }
    }
};

}

#endif
