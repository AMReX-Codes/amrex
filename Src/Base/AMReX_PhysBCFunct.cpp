
#include <AMReX_PhysBCFunct.H>

namespace amrex {

BndryFunctBase::BndryFunctBase ()
    :
    m_func(0),
    m_func3D(0)
{}

BndryFunctBase::BndryFunctBase (BndryFuncDefault inFunc)
    :
    m_func(inFunc),
    m_func3D(0)
{}

BndryFunctBase::BndryFunctBase (BndryFunc3DDefault inFunc)
    :
    m_func(0),
    m_func3D(inFunc)
{}

BndryFunctBase::~BndryFunctBase () {}

BndryFunctBase*
BndryFunctBase::clone () const
{
    return new BndryFunctBase(*this);
}

void
BndryFunctBase::operator () (Real* data,const int* lo,const int* hi,
			     const int* dom_lo, const int* dom_hi,
			     const Real* dx, const Real* grd_lo,
			     const Real* time, const int* bc) const
{
    BL_ASSERT(m_func != 0 || m_func3D != 0);

    if (m_func != 0) {
	m_func(data,ARLIM(lo),ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
    } else {
	m_func3D(data,ARLIM_3D(lo),ARLIM_3D(hi),ARLIM_3D(dom_lo),ARLIM_3D(dom_hi),
		 ZFILL(dx),ZFILL(grd_lo),time,bc);
    }
}


PhysBCFunct::PhysBCFunct (const Geometry& geom, const BCRec& bcr, const BndryFunctBase& func)
    : m_geom(geom), m_bcr(bcr), m_bc_func(func.clone())
{ }

void
PhysBCFunct::define (const Geometry& geom, const BCRec& bcr, const BndryFunctBase& func)
{
    m_geom = geom;
    m_bcr = bcr;
    m_bc_func.reset(func.clone());
}

void
PhysBCFunct::FillBoundary (MultiFab& mf, int, int, Real time)
{
    BL_PROFILE("PhysBCFunct::FillBoundary");

    if (mf.nGrow() == 0) return;
    
    if (m_geom.isAllPeriodic()) return;

    const Box&     domain      = m_geom.Domain();
    const int*     dlo         = domain.loVect();
    const int*     dhi         = domain.hiVect();
    const Real*    dx          = m_geom.CellSize();
    const RealBox& prob_domain = m_geom.ProbDomain();
    const Real*    problo      = prob_domain.lo();

    // create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	if (m_geom.isPeriodic(i)) {
	    gdomain.grow(i, mf.nGrow());
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
	FArrayBox& dest = mf[mfi];
	const Box& bx = dest.box();
        
        // if there are cells not in the valid + periodic grown box
        // we need to fill them here
	if (!gdomain.contains(bx)) {

	    const int* fablo = bx.loVect();
	    const int* fabhi = bx.hiVect();

	    Real xlo[AMREX_SPACEDIM];
	    for (int i = 0; i < AMREX_SPACEDIM; i++)
	    {
		xlo[i] = problo[i] + dx[i]*(fablo[i]-dlo[i]);
	    }

	    (*m_bc_func)(dest.dataPtr(), fablo, fabhi, dlo, dhi,
			 dx, xlo, &time, m_bcr.vect());
	}
    }
}

}
