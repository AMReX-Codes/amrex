
#include "interpolator.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FACINT2  acint2_
#  define FANINT2  anint2_
#else
#  define FACINT2  ACINT2
#  define FANINT2  ANINT2
#endif

extern "C" {
  void FACINT2(Real*, intS, intS, Real*, intS, intS, intRS);
  void FANINT2(Real*, intS, intS, Real*, intS, intS, intRS);
}

Box bilinear_interpolator_class::box(const Box& region,
				     const IntVect& rat) const
{
  if (region.cellCentered()) {
    return grow(coarsen(region, rat), 1);
  }
  else if (region.type() == nodevect) {
    return coarsen(region, rat);
  }
  else {
    BoxLib::Error("bilinear_interpolator_class::box---Interpolation only defined for pure CELL- or NODE-based data");
    return Box();
  }
}

void bilinear_interpolator_class::fill(Fab& patch,
				       const Box& region,
				       Fab& cgr,
				       const Box& cb,
				       const IntVect& rat) const
{
  if (patch.box().cellCentered()) {
    for (int i = 0; i < patch.nComp(); i++) {
      FACINT2(patch.dataPtr(i), dimlist(patch.box()), dimlist(region),
	      cgr.dataPtr(i), dimlist(cgr.box()), dimlist(cb),
	      D_DECL(rat[0], rat[1], rat[2]));
    }
  }
  else if (patch.box().type() == nodevect) {
    Box eregion = refine(cb, rat);
    if (eregion == region) {
      for (int i = 0; i < patch.nComp(); i++) {
	FANINT2(patch.dataPtr(i), dimlist(patch.box()), dimlist(region),
		cgr.dataPtr(i), dimlist(cgr.box()), dimlist(cb),
		D_DECL(rat[0], rat[1], rat[2]));
      }
    }
    else {
      Fab epatch(eregion, patch.nComp());
      for (int i = 0; i < patch.nComp(); i++) {
	FANINT2(epatch.dataPtr(i), dimlist(epatch.box()), dimlist(eregion),
		cgr.dataPtr(i), dimlist(cgr.box()), dimlist(cb),
		D_DECL(rat[0], rat[1], rat[2]));
      }
      patch.copy(epatch,region);
    }
  }
  else
    BoxLib::Error("bilinear_interpolator_class::fill---Interpolation only defined for pure CELL- or NODE-based data");
}
