#ifndef _AMR_DEFS_H_
#define _AMR_DEFS_H_

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

extern bool HG_is_debugging;

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

inline Real mfnorm_0(const MultiFab& mf) 
{
    Real r = 0;
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi ) 
    {
	Real s = cmfi->norm(cmfi->box(), 0, 0, cmfi->nComp());
	r = (r > s) ? r : s;
    }
    ParallelDescriptor::ReduceRealMax(r);
    return r;
}

inline Real mfnorm_2(const MultiFab& mf) 
{
    Real r = 0;
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi ) 
    {
	Real s = cmfi->norm(cmfi->box(), 2, 0, cmfi->nComp());
	r += s*s;
    }
    ParallelDescriptor::ReduceRealSum(r);
    return ::sqrt(r);
}

inline Real mfnorm(const MultiFab& mf) 
{
    Real r = 0;
    assert(mf.nComp() == 1);
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi ) 
    {
	Real s = cmfi->norm(0);
	r = (r > s) ? r : s;
    }
    ParallelDescriptor::ReduceRealMax(r);
    return r;
}

void hg_debug_norm_2(const MultiFab& d, const char* str1, const char* str2);

#ifdef HG_DEBUG
#include <fstream>
extern std::ofstream debug_out;

#define HG_DEBUG_OUT(s) debug_out << s
#define HG_TEST_NORM(d, s) hg_debug_norm_2(d, s, #d)
#else
#define HG_DEBUG_OUT(s)
#define HG_TEST_NORM(d, s)
#endif

#endif
