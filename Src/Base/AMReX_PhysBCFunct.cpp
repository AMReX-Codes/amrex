
#include <AMReX_PhysBCFunct.H>
#include <AMReX_filcc_f.H>

namespace amrex {

void
BndryFuncArray::operator () (Box const& /*bx*/, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp, // BCRec for this box
                             const int /*orig_comp*/)
{
    BL_ASSERT(m_func != nullptr || m_func3D != nullptr);

    const int*  lo = dest.loVect();
    const int*  hi = dest.hiVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    const int* dom_hi = domain.hiVect();
    const Real* dx = geom.CellSize();

    Real grd_lo[AMREX_SPACEDIM];
    const Real* problo = geom.ProbLo();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        grd_lo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }

    static_assert(sizeof(BCRec) == 2*AMREX_SPACEDIM*sizeof(int),
                  "Let us know if this assertion fails");

    for (int icomp = 0; icomp < numcomp; ++icomp)
    {
        Real* data = dest.dataPtr(dcomp+icomp);
        if (m_func != nullptr) {
            m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),
                   dom_lo,dom_hi,
                   dx,grd_lo,&time,bcr[bcomp+icomp].vect());
        } else {
            m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
                     AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                     AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcr[bcomp+icomp].vect());
        }
    }
}

void
CpuBndryFuncFab::operator() (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp,
                             const int orig_comp)
{
    const int* lo = dest.loVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    const Real* dx = geom.CellSize();
    const Real* problo = geom.ProbLo();
    Real xlo[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }
    amrex_fab_filcc(BL_TO_FORTRAN_N_ANYD(dest,dcomp), &numcomp,
                    BL_TO_FORTRAN_BOX(domain),
                    dx, xlo, bcr[bcomp].vect());

    if (f_user != nullptr)
    {
        f_user(bx, dest.array(), dcomp, numcomp, geom.data(), time,
               &(bcr[bcomp]), 0, orig_comp);
    }
}

}
