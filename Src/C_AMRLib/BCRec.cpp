//BL_COPYRIGHT_NOTICE

//
// $Id: BCRec.cpp,v 1.4 1997-12-11 23:27:48 lijewski Exp $
//

#include <BCRec.H>

// -------------------------------------------------------------
BCRec::BCRec(D_DECL(int loX, int loY, int loZ),
             D_DECL(int hiX, int hiY, int hiZ))
{
    D_EXPR(bc[0] = loX,  bc[1] = loY,  bc[2] = loZ);
    D_EXPR(bc[BL_SPACEDIM] = hiX,  bc[BL_SPACEDIM+1] = hiY,  bc[BL_SPACEDIM+2] = hiZ);
}

// -------------------------------------------------------------
BCRec::BCRec(const int* lo, const int* hi)
{
   int i;
   for (i = 0; i < BL_SPACEDIM; i++) {
      bc[i] = lo[i];
      bc[i+BL_SPACEDIM] = hi[i];
   }
}

// -------------------------------------------------------------
BCRec::BCRec(const Box& bx, const Box& domain, const BCRec &bc_domain) 
{
   const int* bxlo = bx.loVect();
   const int* bxhi = bx.hiVect();
   const int* dlo = domain.loVect();
   const int* dhi = domain.hiVect();
   int dir;
   for (dir = 0; dir < BL_SPACEDIM; dir++) {
      int lo = dir;
      int hi = dir+BL_SPACEDIM;
      bc[lo] = ( bxlo[dir]<=dlo[dir] ? bc_domain.bc[lo] : INT_DIR );
      bc[hi] = ( bxhi[dir]>=dhi[dir] ? bc_domain.bc[hi] : INT_DIR );
   }
}

// -------------------------------------------------------------
void
setBC(const Box& bx, const Box& domain, int src_comp, int dest_comp,
      int ncomp, const Array<BCRec> &bc_dom, Array<BCRec> &bcr)
{
   const int* bxlo = bx.loVect();
   const int* bxhi = bx.hiVect();
   const int* dlo = domain.loVect();
   const int* dhi = domain.hiVect();
   int i;
   for (i = 0; i < ncomp; i++) {
      int dc = dest_comp + i;
      int sc = src_comp + i;
      int dir;
      for (dir = 0; dir < BL_SPACEDIM; dir++) {
         bcr[dc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                            ? bc_dom[sc].lo(dir) : INT_DIR ));
         bcr[dc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                            ? bc_dom[sc].hi(dir) : INT_DIR ));
      }
   }
}           

// -------------------------------------------------------------
void
setBC(const Box& bx, const Box& domain, 
      const BCRec &bc_dom, BCRec &bcr)
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo = domain.loVect();
    const int* dhi = domain.hiVect();
    int dir;
    for (dir = 0; dir < BL_SPACEDIM; dir++) {
        bcr.setLo(dir, ( bxlo[dir]<=dlo[dir] ? bc_dom.lo(dir) : INT_DIR ));
        bcr.setHi(dir, ( bxhi[dir]>=dhi[dir] ? bc_dom.hi(dir) : INT_DIR ));
    }
}           

ostream&
operator<< (ostream& os, const BCRec& b)
{
    os << "(BCREC ";
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
        os << b.bc[i] << ":" << b.bc[i+BL_SPACEDIM] << ' ';
    }
    os << ')' << flush;
    return os;
}
