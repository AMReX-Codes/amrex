#include <AMReX_FilCC_C.H>

namespace amrex {

AMREX_GPU_DEVICE
void
filcc_cell (const IntVect& iv, FArrayBox& dest_fab,
            const int dcomp, const int numcomp,
            GeometryData const& geom, const Real time,
            const BCRec* bcr, const int bcomp,
            const int orig_comp)
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

    // xxxxx TODO This does NOT properly fill corner cells.
    
    for (int n = 0; n < numcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n];

        if (i < ilo)
        {
            if (bc.lo(0) == BCType::foextrap)
            {
                q(0,0,0,n) = q(ilo-i,0,0,n);
            }
            else if (bc.lo(0) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.lo(0) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(2*(ilo-i)-1,0,0,n);
            }
            else if (bc.lo(0) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(2*(ilo-i)-1,0,0,n);
            }
        }
        else if (i > ihi)
        {
            if (bc.hi(0) == BCType::foextrap)
            {
                q(0,0,0,n) = q(ihi-i,0,0,n);
            }
            else if (bc.hi(0) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.hi(0) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(2*(ihi-i)+1,0,0,n);
            }
            else if (bc.hi(0) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(2*(ihi-i)+1,0,0,n);
            }
        }

        if (j < jlo)
        {
            if (bc.lo(1) == BCType::foextrap)
            {
                q(0,0,0,n) = q(0,jlo-j,0,n);
            }
            else if (bc.lo(1) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.lo(1) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(0,2*(jlo-j)-1,0,n);
            }
            else if (bc.lo(1) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(0,2*(jlo-j)-1,0,n);
            }
        }
        else if (j > jhi)
        {
            if (bc.hi(1) == BCType::foextrap)
            {
                q(0,0,0,n) = q(0,jhi-j,0,n);
            }
            else if (bc.hi(1) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.hi(1) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(0,2*(jhi-j)+1,0,n);
            }
            else if (bc.hi(1) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(0,2*(jhi-j)+1,0,n);
            }
        }

        if (k < klo)
        {
            if (bc.lo(2) == BCType::foextrap)
            {
                q(0,0,0,n) = q(0,0,klo-k,n);
            }
            else if (bc.lo(2) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.lo(2) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(0,0,2*(klo-k)-1,n);
            }
            else if (bc.lo(2) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(0,0,2*(klo-k)-1,n);
            }
        }
        else if (k > khi)
        {
            if (bc.hi(2) == BCType::foextrap)
            {
                q(0,0,0,n) = q(0,0,khi-k,n);
            }
            else if (bc.hi(2) == BCType::hoextrap)
            {
                amrex::Abort("BCType::hoextrap TODO");
            }
            else if (bc.hi(2) == BCType::reflect_even)
            {
                q(0,0,0,n) = q(0,0,2*(khi-k)+1,n);
            }
            else if (bc.hi(2) == BCType::reflect_odd)
            {
                q(0,0,0,n) = -q(0,0,2*(khi-k)+1,n);
            }
        }

        // xxxxx TODO more hoextrap stuff
    }
}

}
