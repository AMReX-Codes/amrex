
#include <iostream>

#include <AMReX_BCRec.H>

namespace amrex {

void
setBC (const Box&           bx,
       const Box&           domain,
       int                  src_comp,
       int                  dest_comp,
       int                  ncomp,
       const Vector<BCRec>& bc_dom,
       Vector<BCRec>&       bcr)
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int i = 0; i < ncomp; i++)
    {
        int dc = dest_comp + i;
        int sc =  src_comp + i;
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            bcr[dc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                                 ? bc_dom[sc].lo(dir) : BCType::int_dir ));
            bcr[dc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                                 ? bc_dom[sc].hi(dir) : BCType::int_dir ));
        }
    }
}           

void
setBC (const Box&   bx,
       const Box&   domain, 
       const BCRec& bc_dom,
       BCRec&       bcr)
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        bcr.setLo(dir, ( bxlo[dir]<=dlo[dir] ? bc_dom.lo(dir) : BCType::int_dir ));
        bcr.setHi(dir, ( bxhi[dir]>=dhi[dir] ? bc_dom.hi(dir) : BCType::int_dir ));
    }
}           

std::ostream&
operator<< (std::ostream& os,
            const BCRec&  b)
{
    os << "(BCREC ";
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        os << b.bc[i] << ':' << b.bc[i+AMREX_SPACEDIM] << ' ';
    }
    os << ')';
    return os;
}

}
