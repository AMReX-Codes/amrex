//BL_COPYRIGHT_NOTICE

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

#include "hgparallel.H"

extern bool HG_is_debugging;

//
// Index type helper functions:
//

inline IntVect type (const FArrayBox& f) 
{
    return f.box().type();
}

inline IndexType::CellIndex type (const FArrayBox& f, int idir) 
{
    return f.box().type(idir);
}

inline IntVect type (const FabArray<Real,FArrayBox>& r) 
{
    BL_ASSERT(r.length() > 0 && r.boxArray().ok());
    return r.box(0).type();
}

inline IndexType::CellIndex type (const FabArray<Real,FArrayBox>& r, int idir) 
{
    BL_ASSERT(r.length() > 0 && r.boxArray().ok());
    return r.box(0).type(idir);
}

#include "hgdebug_F.H"

#ifdef HG_DEBUG
void hg_debug_norm_2 (const MultiFab& d, const char* str1, const char* str2);

#include <typeinfo>
#include <fstream>
extern std::ofstream debug_out;

#define HG_DEBUG_OUT(s) debug_out << s
#define HG_TEST_NORM(d, s) hg_debug_norm_2(d, s, #d)
#else
#define HG_DEBUG_OUT(s)
#define HG_TEST_NORM(d, s)
#endif

#endif /*_AMR_DEFS_H_*/
