//BL_COPYRIGHT_NOTICE

#include "interpolator.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define FORT_FACINT2  acint2_
#define FORT_FANINT2  anint2_
#elif defined( BL_FORT_USE_UPPERCASE )
#define FORT_FACINT2  ACINT2
#define FORT_FANINT2  ANINT2
#elif defined( BL_FORT_USE_LOWERCASE )
#define FORT_FACINT2  acint2
#define FORT_FANINT2  anint2
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_FACINT2(Real*, intS, intS, const Real*,
		      intS, intS, intRS, const int&);
    void FORT_FANINT2(Real*, intS, intS, const Real*,
		      intS, intS, intRS, const int&);
}

amr_interpolator::~amr_interpolator () {}

Box
bilinear_interpolator::box (const Box&     region,
                                  const IntVect& rat) const
{
    if (region.cellCentered())
    {
	return ::grow(coarsen(region, rat), 1);
    }
    else if (region.type() == IntVect::TheNodeVector())
    {
	return ::coarsen(region, rat);
    }
    else
    {
	BoxLib::Abort( "bilinear_interpolator::box():"
		       "Interpolation only defined for pure CELL- or NODE-based data" ); /*NOTREACHED*/
	return Box();
    }
}

void
bilinear_interpolator::fill (FArrayBox&       patch,
			     const Box&       region,
			     const FArrayBox& cgr,
			     const Box&       cb,
			     const IntVect&   rat) const
{
    if (patch.box().cellCentered())
    {
	FORT_FACINT2(
	    patch.dataPtr(), DIMLIST(patch.box()),
	    DIMLIST(region),
	    cgr.dataPtr(), DIMLIST(cgr.box()),
	    DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
    }
    else if (patch.box().type() == IntVect::TheNodeVector())
    {
        Box eregion = ::refine(cb, rat);

	if (eregion == region)
	{
	    FORT_FANINT2(
		patch.dataPtr(), DIMLIST(patch.box()),
		DIMLIST(region),
		cgr.dataPtr(), DIMLIST(cgr.box()),
		DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
	}
	else
	{
	    FArrayBox epatch(eregion, patch.nComp());
	    FORT_FANINT2(
		epatch.dataPtr(), DIMLIST(epatch.box()),
		DIMLIST(eregion),
		cgr.dataPtr(), DIMLIST(cgr.box()),
		DIMLIST(cb), D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
	    patch.copy(epatch, region);
	}
    }
    else
    {
	BoxLib::Abort( "bilinear_interpolator::fill():"
		       "Interpolation only defined for pure CELL- or NODE-based data" );
    }
}
