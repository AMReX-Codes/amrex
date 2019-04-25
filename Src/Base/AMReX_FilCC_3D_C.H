#ifndef AMREX_FILCC_3D_C_H_
#define AMREX_FILCC_3D_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

AMREX_GPU_DEVICE
inline
void
filcc_cell (const IntVect& iv, FArrayBox& dest_fab,
            const int dcomp, const int numcomp,
            GeometryData const& geom, const Real time,
            const BCRec* bcr, const int bcomp,
            const int orig_comp) noexcept
{
    const int i = iv[0];
    const int j = iv[1];
    const int k = iv[2];
    const auto q = dest_fab.view(iv,dcomp);

    const Box& domain_box = geom.Domain();
    const auto& domain_lo = domain_box.loVect();
    const auto& domain_hi = domain_box.hiVect();
    const int ilo = domain_lo[0];
    const int jlo = domain_lo[1];
    const int klo = domain_lo[2];
    const int ihi = domain_hi[0];
    const int jhi = domain_hi[1];
    const int khi = domain_hi[2];
    const int* qlo = dest_fab.loVect();
    const int* qhi = dest_fab.hiVect();
    const int is = amrex::max(qlo[0],ilo);
    const int js = amrex::max(qlo[1],jlo);
    const int ks = amrex::max(qlo[2],klo);
    const int ie = amrex::min(qhi[0],ihi);
    const int je = amrex::min(qhi[1],jhi);
    const int ke = amrex::min(qhi[2],khi);

    for (int n = 0; n < numcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n];

        if (i < ilo)
        {
            switch (bc.lo(0)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(ilo-i,0,0,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (i < ilo - 1)
                {
                    q(0,0,0,n) = q(ilo-i,0,0,n);
                }
                // i == ilo-1
                else if (ilo+2 <= ie)
                {
                    q(0,0,0,n) = 0.125*(15.*q(1,0,0,n) - 10.*q(2,0,0,n) + 3.*q(3,0,0,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(1,0,0,n) - q(2,0,0,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(2*(ilo-i)-1,0,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(2*(ilo-i)-1,0,0,n);
                break;
            }
            }
        }
        else if (i > ihi)
        {
            switch (bc.hi(0)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(ihi-i,0,0,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (i > ihi + 1)
                {
                    q(0,0,0,n) = q(ihi-i,0,0,n);
                }
                // i == ihi+1
                else if (ihi-2 >= is)
                {
                    q(0,0,0,n) = 0.125*(15.*q(-1,0,0,n) - 10.*q(-2,0,0,n) + 3.*q(-3,0,0,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(-1,0,0,n) - q(-2,0,0,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(2*(ihi-i)+1,0,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(2*(ihi-i)+1,0,0,n);
                break;
            }
            }
        }

        if (j < jlo)
        {
            switch (bc.lo(1)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(0,jlo-j,0,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (j < jlo - 1)
                {
                    q(0,0,0,n) = q(0,jlo-j,0,n);
                }
                // j == jlo-1
                else if (jlo+2 <= je)
                {
                    q(0,0,0,n) = 0.125*(15.*q(0,1,0,n) - 10.*q(0,2,0,n) + 3.*q(0,3,0,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(0,1,0,n) - q(0,2,0,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(0,2*(jlo-j)-1,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(0,2*(jlo-j)-1,0,n);
                break;
            }
            }
        }
        else if (j > jhi)
        {
            switch (bc.hi(1)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(0,jhi-j,0,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (j > jhi + 1)
                {
                    q(0,0,0,n) = q(0,jhi-j,0,n);
                }
                // j == jhi+1
                else if (jhi-2 >= js)
                {
                    q(0,0,0,n) = 0.125*(15.*q(0,-1,0,n) - 10.*q(0,-2,0,n) + 3.*q(0,-3,0,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(0,-1,0,n) - q(0,-2,0,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(0,2*(jhi-j)+1,0,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(0,2*(jhi-j)+1,0,n);
                break;
            }
            }
        }

        if (k < klo)
        {
            switch (bc.lo(2)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(0,0,klo-k,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (k < klo - 1)
                {
                    q(0,0,0,n) = q(0,0,klo-k,n);
                }
                // k == klo-1
                else if (klo+2 <= ke)
                {
                    q(0,0,0,n) = 0.125*(15.*q(0,0,1,n) - 10.*q(0,0,2,n) + 3.*q(0,0,3,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(0,0,1,n) - q(0,0,2,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(0,0,2*(klo-k)-1,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(0,0,2*(klo-k)-1,n);
                break;
            }
            }
        }
        else if (k > khi)
        {
            switch (bc.hi(2)) {
            case (BCType::foextrap):
            {
                q(0,0,0,n) = q(0,0,khi-k,n);
                break;
            }
            case (BCType::hoextrap):
            {
                if (k > khi + 1)
                {
                    q(0,0,0,n) = q(0,0,khi-k,n);
                }
                // k == khi+1
                else if (khi-2 >= ks)
                {
                    q(0,0,0,n) = 0.125*(15.*q(0,0,-1,n) - 10.*q(0,0,-2,n) + 3.*q(0,0,-3,n));
                }
                else
                {
                    q(0,0,0,n) = 0.5*(3.*q(0,0,-1,n) - q(0,0,-2,n));
                }
                break;
            }
            case (BCType::reflect_even):
            {
                q(0,0,0,n) = q(0,0,2*(khi-k)+1,n);
                break;
            }
            case (BCType::reflect_odd):
            {
                q(0,0,0,n) = -q(0,0,2*(khi-k)+1,n);
                break;
            }
            }
        }
    }
}

}

#endif
