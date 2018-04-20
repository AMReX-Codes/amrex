
#include <iostream>

#include <AMReX_BCRec.H>

namespace amrex {

BCRec::BCRec ()
    : bc {AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus),
          AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}
{ }

BCRec::BCRec (AMREX_D_DECL(int loX, int loY, int loZ),
              AMREX_D_DECL(int hiX, int hiY, int hiZ))
{
    AMREX_D_EXPR(bc[0] = loX,  bc[1] = loY,  bc[2] = loZ);
    AMREX_D_EXPR(bc[AMREX_SPACEDIM]=hiX,  bc[AMREX_SPACEDIM+1]=hiY,  bc[AMREX_SPACEDIM+2]=hiZ);
}

BCRec::BCRec (const int* a_lo,
              const int* a_hi)
{
    BL_ASSERT(!(a_lo == 0));
    BL_ASSERT(!(a_hi == 0));

    AMREX_D_TERM(bc[0] = a_lo[0];,
                 bc[1] = a_lo[1];,
                 bc[2] = a_lo[2];);

    AMREX_D_TERM(bc[AMREX_SPACEDIM+0] = a_hi[0];,
                 bc[AMREX_SPACEDIM+1] = a_hi[1];,
                 bc[AMREX_SPACEDIM+2] = a_hi[2];);
}

BCRec::BCRec (const Box&   bx,
              const Box&   domain,
              const BCRec& bc_domain) 
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        int ilo = dir;
        int ihi = dir+AMREX_SPACEDIM;
        bc[ilo] = ( bxlo[dir]<=dlo[dir] ? bc_domain.bc[ilo] : BCType::int_dir );
        bc[ihi] = ( bxhi[dir]>=dhi[dir] ? bc_domain.bc[ihi] : BCType::int_dir );
    }
}

void
setBC (const Box&          bx,
               const Box&          domain,
               int                 src_comp,
               int                 dest_comp,
               int                 ncomp,
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
