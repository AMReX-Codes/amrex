#ifndef _AMR_DEFS_H_
#define _AMR_DEFS_H_

// includes IntVect.H, IndexType.H, Real.H, FArrayBox.H, BoxArray.H etc.
#include <MultiFab.H>

#undef SPACEDIM

#if (BL_SPACEDIM == 1)
#  define intRS const int&
#  define intS  const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0)
#  define RealRS const Real&
#  define RealPS Real*
#elif (BL_SPACEDIM == 2)
#  define intRS const int&, const int&
#  define intS  const int&, const int&, const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1)
#  define RealRS const Real&, const Real&
#  define RealPS Real*, Real*
#elif (BL_SPACEDIM == 3)
#  define intRS const int&, const int&, const int&
#  define intS  const int&, const int&, const int&, const int&, const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1), a.smallEnd(2), a.bigEnd(2)
#  define RealRS const Real&, const Real&, const Real&
#  define RealPS Real*, Real*, Real*
#endif

// index type helper functions:

inline IntVect type(const FArrayBox& f) {
  return f.box().type();
}

inline IndexType::CellIndex type(const FArrayBox& f, int idir) {
  return f.box().type(idir);
}

inline IntVect type(const FabArray<Real,FArrayBox>& r) {
  assert(r.length() > 0 && r.boxArray().ok());
  return r.box(0).type();
}

inline IndexType::CellIndex type(const FabArray<Real,FArrayBox>& r, int idir) {
  assert(r.length() > 0 && r.boxArray().ok());
  return r.box(0).type(idir);
}

// norm helper functions:

inline Real mfnorm(MultiFab& mf) {
  Real r = 0;
  int i;
  for (i = 0; i < mf.length(); i++) {
    Real s = mf[i].norm(0);
    r = (r > s) ? r : s;
  }
  return r;
}

inline Real pmfnorm(PArray<MultiFab>& pmf) {
  Real r = 0;
  int i;
  for (i = 0; i < pmf.length(); i++) {
    Real s = mfnorm(pmf[i]);
    r = (r > s) ? r : s;
  }
  return r;
}

inline Real pmfmin(PArray<MultiFab>& pmf) {
  Real r = pmf[0].min(0);
  int i;
  for (i = 1; i < pmf.length(); i++) {
    Real s = pmf[i].min(0);
    r = (r < s) ? r : s;
  }
  return r;
}

inline Real pmfmax(PArray<MultiFab>& pmf) {
  Real r = pmf[0].max(0);
  int i;
  for (i = 1; i < pmf.length(); i++) {
    Real s = pmf[i].max(0);
    r = (r > s) ? r : s;
  }
  return r;
}

#endif
