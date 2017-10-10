
#include <AMReX_MultiMask.H>
#include <AMReX_BndryData.H>

namespace amrex {

MultiMask::MultiMask (const BoxArray& ba, const DistributionMapping& dm, int ncomp)
    : m_fa(ba, dm, ncomp, 0, MFInfo(), DefaultFabFactory<Mask>())
{ }

MultiMask::MultiMask (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
		      Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval)
{
    define(regba, dm, geom, face, in_rad, out_rad, extent_rad, ncomp, initval);
}

void
MultiMask::define (const BoxArray& ba, const DistributionMapping& dm, int ncomp)
{
    BL_ASSERT(m_fa.size() == 0);
    m_fa.define(ba,dm,ncomp,0,MFInfo(),DefaultFabFactory<Mask>());
}

void
MultiMask::define (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
		   Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval)
{
    BL_ASSERT(m_fa.size() == 0);

    BndryBATransformer bbatrans(face,IndexType::TheCellType(),in_rad,out_rad,extent_rad);
    BoxArray mskba(regba, bbatrans);
    m_fa.define(mskba, dm, ncomp, 0, MFInfo(), DefaultFabFactory<Mask>());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    if (initval)
    {
	Vector<IntVect> pshifts(26);
	std::vector< std::pair<int,Box> > isects;

	for (MFIter mfi(m_fa); mfi.isValid(); ++mfi)
	{
	    Mask& m = m_fa[mfi];
	    const Box& face_box = m.box();
		
	    m.setVal(BndryData::outside_domain);
	    const Box& dbox = geom.Domain() & face_box;
	    m.setVal(BndryData::not_covered,dbox,0,ncomp);
	    //
	    // Now have to set as not_covered the periodic translates as well.
	    //
	    if (geom.isAnyPeriodic() && !geom.Domain().contains(face_box))
	    {
		geom.periodicShift(geom.Domain(), face_box, pshifts);
		
		for (Vector<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
		     it != End;
		     ++it)
		{
		    const IntVect& iv = *it;
		    m.shift(iv);
		    const Box& target = geom.Domain() & m.box();
		    m.setVal(BndryData::not_covered,target,0,ncomp);
		    m.shift(-iv);
		}
	    }
	    //
	    // Turn mask off on intersection with regba
	    //
	    regba.intersections(face_box,isects);
	    
	    for (int ii = 0, N = isects.size(); ii < N; ii++) {
		m.setVal(BndryData::covered, isects[ii].second, 0, ncomp);
	    }
	    
	    if (geom.isAnyPeriodic() && !geom.Domain().contains(face_box))
	    {
		//
		// Handle special cases if periodic: "face_box" hasn't changed;
		// reuse pshifts from above.
		//
		for (Vector<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
		     it != End;
		     ++it)
		{
		    const IntVect& iv = *it;
		    m.shift(iv);
		    regba.intersections(m.box(),isects);
		    for (int ii = 0, N = isects.size(); ii < N; ii++) {
			m.setVal(BndryData::covered, isects[ii].second, 0, ncomp);
		    }
		    m.shift(-iv);
		}
	    }
	}
    }
}

void 
MultiMask::Copy (MultiMask& dst, const MultiMask& src)
{
    BL_ASSERT(dst.nComp() == src.nComp());
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.DistributionMap() == src.DistributionMap());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst.m_fa); mfi.isValid(); ++mfi) {
	dst.m_fa[mfi].copy(src.m_fa[mfi]);
    }
}

}
