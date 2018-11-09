
#include <AMReX_PhysBCFunct.H>

namespace amrex {

void
BndryFuncArray::operator () (Box const& /*bx*/, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp, // BCRec for this box
                             const int /*orig_comp*/)
{
    BL_ASSERT(m_func != nullptr || m_func3D != nullptr);

    Real*     data = dest.dataPtr(dcomp);
    const Box&  bx = dest.box();
    const int*  lo = dest.loVect();
    const int*  hi = dest.hiVect();
    const Box& domain = geom.Domain();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();
    const Real* dx = geom.CellSize();

    Real grd_lo[AMREX_SPACEDIM];
    const Real* problo = geom.ProbLo();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        grd_lo[i] = problo[i] + dx[i]*(lo[i]-plo[i]);
    }

    static_assert(sizeof(BCRec) == 2*AMREX_SPACEDIM*sizeof(int),
                  "Let us know if this assertion fails");

    if (m_func != nullptr) {
	m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),plo,phi,dx,grd_lo,&time,bcr[0].vect());
    } else {
	m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(plo),AMREX_ARLIM_3D(phi),
		 AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcr[0].vect());
    }
}

}
