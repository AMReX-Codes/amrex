#ifndef AMREX_FILCC_3D_C_H_
#define AMREX_FILCC_3D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

struct FilccCell
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

            if (i < ilo)
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
                    else if (ilo+2 <= amrex::min(q.end.x-1,ihi))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i+1,j,k,n) - Real(10.)*q(i+2,j,k,n) + Real(3.)*q(i+3,j,k,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i+1,j,k,n) - q(i+2,j,k,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(ilo-i)*(q(ilo,j,k,n) - q(ilo+1,j,k,n)) + q(ilo,j,k,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(2*ilo-i-1,j,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(2*ilo-i-1,j,k,n);
                    break;
                }
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
                        q(i,j,k,n) = q(ihi,j,k,n);
                    }
                    // i == ihi+1
                    else if (ihi-2 >= amrex::max(q.begin.x,ilo))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i-1,j,k,n) - Real(10.)*q(i-2,j,k,n) + Real(3.)*q(i-3,j,k,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i-1,j,k,n) - q(i-2,j,k,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(i-ihi)*(q(ihi,j,k,n) - q(ihi-1,j,k,n)) + q(ihi,j,k,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(2*ihi-i+1,j,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(2*ihi-i+1,j,k,n);
                    break;
                }
                }
            }

            if (j < jlo)
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
                    else if (jlo+2 <= amrex::min(q.end.y-1,jhi))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i,j+1,k,n) - Real(10.)*q(i,j+2,k,n) + Real(3.)*q(i,j+3,k,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i,j+1,k,n) - q(i,j+2,k,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(jlo-j)*(q(i,jlo,k,n) - q(i,jlo+1,k,n)) + q(i,jlo,k,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(i,2*jlo-j-1,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(i,2*jlo-j-1,k,n);
                    break;
                }
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
                    else if (jhi-2 >= amrex::max(q.begin.y,jlo))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i,j-1,k,n) - Real(10.)*q(i,j-2,k,n) + Real(3.)*q(i,j-3,k,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i,j-1,k,n) - q(i,j-2,k,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(j-jhi)*(q(i,jhi,k,n) - q(i,jhi-1,k,n)) + q(i,jhi,k,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(i,2*jhi-j+1,k,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(i,2*jhi-j+1,k,n);
                    break;
                }
                }
            }

            if (k < klo)
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
                    else if (klo+2 <= amrex::min(q.end.z-1,khi))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i,j,k+1,n) - Real(10.)*q(i,j,k+2,n) + Real(3.)*q(i,j,k+3,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i,j,k+1,n) - q(i,j,k+2,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(klo-k)*(q(i,j,klo,n) - q(i,j,klo+1,n)) + q(i,j,klo,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(i,j,2*klo-k-1,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(i,j,2*klo-k-1,n);
                    break;
                }
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
                    else if (khi-2 >= amrex::max(q.begin.z,klo))
                    {
                        q(i,j,k,n) = Real(0.125)*(Real(15.)*q(i,j,k-1,n) - Real(10.)*q(i,j,k-2,n) + Real(3.)*q(i,j,k-3,n));
                    }
                    else
                    {
                        q(i,j,k,n) = Real(0.5)*(Real(3.)*q(i,j,k-1,n) - q(i,j,k-2,n));
                    }
                    break;
                }
                case (BCType::hoextrapcc):
                {
                    q(i,j,k,n) = Real(k-khi)*(q(i,j,khi,n) - q(i,j,khi-1,n)) + q(i,j,khi,n);
                    break;
                }
                case (BCType::reflect_even):
                {
                    q(i,j,k,n) = q(i,j,2*khi-k+1,n);
                    break;
                }
                case (BCType::reflect_odd):
                {
                    q(i,j,k,n) = -q(i,j,2*khi-k+1,n);
                    break;
                }
                }
            }
        }
    }
};

}

#endif
