
#include "interpolator.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FACINT2  acint2_
#define FORT_FANINT2  anint2_
#else
#define FORT_FACINT2  ACINT2
#define FORT_FANINT2  ANINT2
#endif

extern "C" 
{
    void FORT_FACINT2(Real*, intS, intS, const Real*, intS, intS, intRS, const int&);
    void FORT_FANINT2(Real*, intS, intS, const Real*, intS, intS, intRS, const int&);
}

Box bilinear_interpolator_class::box(const Box& region, const IntVect& rat) const
{
    if (region.cellCentered()) 
    {
	return grow(coarsen(region, rat), 1);
    }
    else if (region.type() == IntVect::TheNodeVector()) 
    {
	return coarsen(region, rat);
    }
    else 
    {
	BoxLib::Abort( "bilinear_interpolator_class::box(): Interpolation only defined for pure CELL- or NODE-based data" ); /*NOTREACHED*/
	return Box();
    }
}

void bilinear_interpolator_class::fill(FArrayBox& patch, const Box& region, const FArrayBox& cgr, const Box& cb, const IntVect& rat) const
{
    if (patch.box().cellCentered()) 
    {
	FORT_FACINT2(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
    }
    else if (patch.box().type() == IntVect::TheNodeVector()) 
    {
	const Box eregion = refine(cb, rat);
	if (eregion == region) 
	{
	    FORT_FANINT2(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
	}
	else 
	{
	    FArrayBox epatch(eregion, patch.nComp());
	    FORT_FANINT2(epatch.dataPtr(), DIMLIST(epatch.box()), DIMLIST(eregion), cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
	    patch.copy(epatch,region);
	}
    }
    else
    {
	BoxLib::Abort( "bilinear_interpolator_class::fill(): Interpolation only defined for pure CELL- or NODE-based data" );
    }
}
