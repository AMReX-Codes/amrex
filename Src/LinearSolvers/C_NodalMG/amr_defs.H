#ifndef _AMR_DEFS_H_
#define _AMR_DEFS_H_

//
// $Id: amr_defs.H,v 1.13 1998-05-15 00:31:48 car Exp $
//

#include <MultiFab.H>

#if (BL_SPACEDIM == 1)
#define intRS const int&
#define intS  const int&, const int&
#define DIMLIST(a) a.smallEnd(0), a.bigEnd(0)
#define RealPS Real*
#define CRealPS const Real*
#elif (BL_SPACEDIM == 2)
#define intRS const int&, const int&
#define intS  const int&, const int&, const int&, const int&
#define DIMLIST(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1)
#define RealPS Real*, Real*
#define CRealPS const Real*, const Real*
#elif (BL_SPACEDIM == 3)
#define intRS const int&, const int&, const int&
#define intS  const int&, const int&, const int&, const int&, const int&, const int&
#define DIMLIST(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1), a.smallEnd(2), a.bigEnd(2)
#define RealPS Real*, Real*, Real*
#define CRealPS const Real*, const Real*, const Real*
#endif

#include "hgparallel.h"

// index type helper functions:

inline IntVect type(const FArrayBox& f) 
{
    return f.box().type();
}

inline IndexType::CellIndex type(const FArrayBox& f, int idir) 
{
    return f.box().type(idir);
}

inline IntVect type(const FabArray<Real,FArrayBox>& r) 
{
    assert(r.length() > 0 && r.boxArray().ok());
    return r.box(0).type();
}

inline IndexType::CellIndex type(const FabArray<Real,FArrayBox>& r, int idir) 
{
    assert(r.length() > 0 && r.boxArray().ok());
    return r.box(0).type(idir);
}

// norm helper functions:

inline Real mfnorm(const MultiFab& mf) 
{
    Real r = 0;
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi ) 
    {
	Real s = cmfi->norm(0);
	r = (r > s) ? r : s;
    }
    ParallelDescriptor::ReduceRealMax(r);
    return r;
}

#endif
