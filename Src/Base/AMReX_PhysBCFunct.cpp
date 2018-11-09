
#include <AMReX_PhysBCFunct.H>

namespace amrex {

void
BndryFuncPtr::operator () (Box const& /*bx*/, FArrayBox& dest,
                           const int dcomp, const int numcomp,
                           Geometry const& geom, const Real time,
                           const Vector<BCRec>& bcr, const int bcomp)
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

    BCRec bcrtmp;
    Vector<int> bcrs(2*AMREX_SPACEDIM*numcomp);
    int* bci = bcrs.data();
    for (int j = 0; j < numcomp; ++j)
    {
        amrex::setBC(bx, domain, bcr[bcomp+j], bcrtmp);
        const int* bc = bcrtmp.vect();
        for (int k = 0; k < 2*AMREX_SPACEDIM; k++) {
            bci[k] = bc[k];
        }
        bci += 2*AMREX_SPACEDIM;
    }

    if (m_func != nullptr) {
	m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),plo,phi,dx,grd_lo,&time,bcrs.data());
    } else {
	m_func3d(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(plo),AMREX_ARLIM_3D(phi),
		 AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcrs.data());
    }
}

}
