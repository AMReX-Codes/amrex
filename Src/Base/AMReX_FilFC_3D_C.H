#ifndef AMREX_FILFC_3D_C_H_
#define AMREX_FILFC_3D_C_H_
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
        const int k = iv[2];

        // Domain box is indexed according to the face currently treating
        const IndexType &idxType = domain_box.ixType();
        const auto& domain_lo = domain_box.loVect();
        const auto& domain_hi = domain_box.hiVect();
        const int ilo = domain_lo[0];
        const int jlo = domain_lo[1];
        const int klo = domain_lo[2];
        const int ihi = domain_hi[0];
        const int jhi = domain_hi[1];
        const int khi = domain_hi[2];

        for (int n = dcomp; n < numcomp+dcomp; ++n)
        {
            const BCRec& bc = bcr[bcomp+n-dcomp];

            if (i == ilo)
            {
                // Enforce reflect_odd on the x domain face
                if ((bc.lo(0) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(0)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (i < ilo)
            {
                switch (bc.lo(0)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(ilo,j,k,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i < ilo - 1)
                    {
                        q(i,j,k,n) = q(ilo,j,k,n);
                    }
                    // i == ilo-1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(0)) ? Real(2.0)*q(i+1,j,k,n) - q(i+2,j,k,n)
                                                                : Real(0.5)*(Real(3.)*q(i+1,j,k,n) - q(i+2,j,k,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(0)) ? q(2*ilo-i,j,k,n)
                                                           : q(2*ilo-i-1,j,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(0)) ? -q(2*ilo-i,j,k,n)
                                                           : -q(2*ilo-i-1,j,k,n);
                    break;
                }
                }
            }
            if (i == ihi)
            {
                // Enforce reflect_odd on the x domain face
                if ((bc.hi(0) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(0)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (i > ihi)
            {
                switch (bc.hi(0)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(ihi,j,k,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (i > ihi + 1)
                    {
                        q(i,j,k,n) = q(ilo,j,k,n);
                    }
                    // i == ihi+1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(0)) ? Real(2.0)*q(i-1,j,k,n) - q(i-2,j,k,n)
                                                                : Real(0.5)*(Real(3.)*q(i-1,j,k,n) - q(i-2,j,k,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(0)) ? q(2*ihi-i,j,k,n)
                                                           : q(2*ihi-i+1,j,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(0)) ? -q(2*ihi-i,j,k,n)
                                                           : -q(2*ihi-i+1,j,k,n);
                    break;
                }
                }
            }

            if (j == jlo)
            {
                // Enforce reflect_odd on the y domain face
                if ((bc.lo(1) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(1)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (j < jlo)
            {
                switch (bc.lo(1)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(i,jlo,k,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (j < jlo - 1)
                    {
                        q(i,j,k,n) = q(i,jlo,k,n);
                    }
                    // j == jlo-1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(1)) ? Real(2.0)*q(i,j+1,k,n) - q(i,j+2,k,n)
                                                                : Real(0.5)*(Real(3.)*q(i,j+1,k,n) - q(i,j+2,k,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(1)) ? q(i,2*jlo-j,k,n)
                                                           : q(i,2*jlo-j-1,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(1)) ? -q(i,2*jlo-j,k,n)
                                                           : -q(i,2*jlo-j-1,k,n);
                    break;
                }
                }
            }
            else if (j == jhi)
            {
                // Enforce reflect_odd on the y domain face
                if ((bc.hi(1) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(1)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (j > jhi)
            {
                switch (bc.hi(1)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(i,jhi,k,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (j > jhi + 1)
                    {
                        q(i,j,k,n) = q(i,jhi,k,n);
                    }
                    // j == jhi+1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(1)) ? Real(2.0)*q(i,j-1,k,n) - q(i,j-2,k,n)
                                                                : Real(0.5)*(Real(3.)*q(i,j-1,k,n) - q(i,j-2,k,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(1)) ? q(i,2*jhi-j,k,n)
                                                           : q(i,2*jhi-j+1,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(1)) ? -q(i,2*jhi-j,k,n)
                                                           : -q(i,2*jhi-j+1,k,n);
                    break;
                }
                }
            }

            if (k == klo)
            {
                // Enforce reflect_odd on the z domain face
                if ((bc.lo(2) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(2)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (k < klo)
            {
                switch (bc.lo(2)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(i,j,klo,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (k < klo - 1)
                    {
                        q(i,j,k,n) = q(i,j,klo,n);
                    }
                    // k == klo-1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(2)) ? Real(2.0)*q(i,j,k+1,n) - q(i,j,k+2,n)
                                                                : Real(0.5)*(Real(3.)*q(i,j,k+1,n) - q(i,j,k+2,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(2)) ? q(i,j,2*klo-k,n)
                                                           : q(i,j,2*klo-k-1,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(2)) ? -q(i,j,2*klo-k,n)
                                                           : -q(i,j,2*klo-k-1,n);
                    break;
                }
                }
            }
            if (k == khi)
            {
                // Enforce reflect_odd on the z domain face
                if ((bc.hi(2) == BCType::reflect_odd) &&
                    (idxType.nodeCentered(2)) ) {
                    q(i,j,k,n) = 0.0;
                }
            }
            else if (k > khi)
            {
                switch (bc.hi(2)) {
                case (BCType::foextrap):
                {
                    q(i,j,k,n) = q(i,j,khi,n);
                    break;
                }
                case (BCType::hoextrap):
                {
                    if (k > khi + 1)
                    {
                        q(i,j,k,n) = q(i,j,khi,n);
                    }
                    // k == khi+1
                    else
                    {
                        q(i,j,k,n) =  (idxType.nodeCentered(2)) ? Real(2.0)*q(i,j,k-1,n) - q(i,j,k-2,n)
                                                                : Real(0.5)*(Real(3.)*q(i,j,k-1,n) - q(i,j,k-2,n));
                    }
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(2)) ? q(i,j,2*khi-k,n)
                                                           : q(i,j,2*khi-k+1,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = (idxType.nodeCentered(2)) ? -q(i,j,2*khi-k,n)
                                                           : -q(i,j,2*khi-k+1,n);
                    break;
                }
                }
            }
        }
    }
};

}

#endif
