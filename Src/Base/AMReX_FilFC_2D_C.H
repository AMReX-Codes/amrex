#ifndef AMREX_FILFC_2D_C_H_
#define AMREX_FILFC_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

AMREX_GPU_DEVICE
inline
void
filfc_cell (const IntVect& iv, const IndexType &idxType,
            Array4<Real> const& q,
            const int dcomp, const int numcomp,
            GeometryData const& geom, const Real /*time*/,
            const BCRec* bcr, const int bcomp,
            const int /*orig_comp*/) noexcept
{
    const int i = iv[0];
    const int j = iv[1];

    const Box& domain_box = geom.Domain();
    const auto& domain_lo = domain_box.loVect();
    const auto& domain_hi = domain_box.hiVect();
    const int ilo = domain_lo[0];
    const int jlo = domain_lo[1];
    const int ihi = domain_hi[0];
    const int jhi = domain_hi[1];
    const int is = amrex::max(q.begin.x,ilo);
    const int js = amrex::max(q.begin.y,jlo);
    const int ie = amrex::min(q.end.x-1,ihi);
    const int je = amrex::min(q.end.y-1,jhi);

    for (int n = dcomp; n < numcomp+dcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n-dcomp];

        if (i < ilo)
        {
            switch (bc.lo(0)) {
            case (BCType::foextrap):
            {
                q(i,j,0,n) = q(ilo,j,0,n);
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
        else if (i > ihi)
        {
            switch (bc.hi(0)) {
            case (BCType::foextrap):
            {
                q(i,j,0,n) = q(ihi,j,0,n);
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

        if (j < jlo)
        {
            switch (bc.lo(1)) {
            case (BCType::foextrap):
            {
                q(i,j,0,n) = q(i,jlo,0,n);
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
        else if (j > jhi)
        {
            switch (bc.hi(1)) {
            case (BCType::foextrap):
            {
                q(i,j,0,n) = q(i,jhi,0,n);
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

}

#endif
