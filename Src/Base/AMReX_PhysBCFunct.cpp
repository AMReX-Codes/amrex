
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

    Real*     data = dest.dataPtr(dcomp);
    const Box&  bx = dest.box();
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

    if (m_func != nullptr) {
	m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),
               dom_lo,dom_hi,
               dx,grd_lo,&time,bcr[bcomp].vect());
    } else {
	m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
                 AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
		 AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcr[bcomp].vect());
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
        f_user(BL_TO_FORTRAN_N_ANYD(dest,dcomp), numcomp,
               BL_TO_FORTRAN_BOX(domain),
               dx, problo, time, bcr[bcomp].vect(), orig_comp);
    }
}

void
GpuBndryFuncFab::operator() (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp,
                             const int orig_comp)
{
    FArrayBox* fab = &dest;
    const auto geomdata = geom.data();

    Gpu::AsyncArray<BCRec> bcr_aa(bcr.data()+bcomp, numcomp);
    BCRec* bcr_p = bcr_aa.data();
    
    Box gdomain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            gdomain.grow(idim,bx.length(idim));
        }
    }

    const auto& bl = amrex::boxDiff(bx, gdomain);
    if (bl.isEmpty()) return;

    const Vector<Box>& boxes = bl.data();
    const int nboxes = boxes.size();
    Gpu::AsyncArray<Box> boxes_aa(boxes.data(), nboxes);
    Box* boxes_p = boxes_aa.data();

    long ncells = 0;
    for (const auto& b : boxes) {
        ncells += b.numPts();
    }

    AMREX_LAUNCH_DEVICE_LAMBDA (ncells, icell,
    {
        int ibox;
        long ncells_subtotal = 0;
        long offset = 0;
        for (ibox = 0; ibox < nboxes; ++ibox) {
            const long n = boxes_p[ibox].numPts();
            ncells_subtotal += n;
            if (icell < ncells_subtotal) {
                offset = icell - (ncells_subtotal - n);
                break;
            }
        }

        const Box& b = boxes_p[ibox];
        const IntVect& iv = b.atOffset(offset);

        //        amrex_fill_bc_cc(tbx, *fab, dcomp, numcomp, geomdata, time,
//                         bcr_aa, scomp);
    });
}

}
