
#include <iostream>

#include <BCRec.H>

BCRec::BCRec (D_DECL(int loX, int loY, int loZ),
              D_DECL(int hiX, int hiY, int hiZ))
{
    D_EXPR(bc[0] = loX,  bc[1] = loY,  bc[2] = loZ);
    D_EXPR(bc[BL_SPACEDIM]=hiX,  bc[BL_SPACEDIM+1]=hiY,  bc[BL_SPACEDIM+2]=hiZ);
}

BCRec::BCRec (const int* lo,
              const int* hi)
{
    BL_ASSERT(!(lo == 0));
    BL_ASSERT(!(hi == 0));

    D_TERM(bc[0] = lo[0];,
           bc[1] = lo[1];,
           bc[2] = lo[2];);

    D_TERM(bc[BL_SPACEDIM+0] = hi[0];,
           bc[BL_SPACEDIM+1] = hi[1];,
           bc[BL_SPACEDIM+2] = hi[2];);
}

BCRec::BCRec (const Box&   bx,
              const Box&   domain,
              const BCRec& bc_domain) 
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int lo = dir;
        int hi = dir+BL_SPACEDIM;
        bc[lo] = ( bxlo[dir]<=dlo[dir] ? bc_domain.bc[lo] : INT_DIR );
        bc[hi] = ( bxhi[dir]>=dhi[dir] ? bc_domain.bc[hi] : INT_DIR );
    }
}

void
BoxLib::setBC (const Box&          bx,
               const Box&          domain,
               int                 src_comp,
               int                 dest_comp,
               int                 ncomp,
               const Array<BCRec>& bc_dom,
               Array<BCRec>&       bcr)
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int i = 0; i < ncomp; i++)
    {
        int dc = dest_comp + i;
        int sc =  src_comp + i;
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            bcr[dc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                                 ? bc_dom[sc].lo(dir) : INT_DIR ));
            bcr[dc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                                 ? bc_dom[sc].hi(dir) : INT_DIR ));
        }
    }
}           

void
BoxLib::setBC (const Box&   bx,
               const Box&   domain, 
               const BCRec& bc_dom,
               BCRec&       bcr)
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        bcr.setLo(dir, ( bxlo[dir]<=dlo[dir] ? bc_dom.lo(dir) : INT_DIR ));
        bcr.setHi(dir, ( bxhi[dir]>=dhi[dir] ? bc_dom.hi(dir) : INT_DIR ));
    }
}           

std::ostream&
operator<< (std::ostream& os,
            const BCRec&  b)
{
    os << "(BCREC ";
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        os << b.bc[i] << ':' << b.bc[i+BL_SPACEDIM] << ' ';
    }
    os << ')';
    return os;
}
